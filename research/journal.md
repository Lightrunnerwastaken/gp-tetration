# Experiment-Journal: fatou_fork.gp

Format pro Eintrag:

## exp-NNN (YYYY-MM-DD) — Kurztitel
- **Hypothese:**
- **Mutation:** (Datei + was)
- **Gate:** PASS/FAIL
- **Benchmark:** Label, warm-Median digits/s vs. Vorgaenger, Ergebnisdatei
- **Entscheidung:** keep/revert
- **Learnings:**

---

## exp-000 (2026-07-10) — Baseline
- **Hypothese:** — (Ausgangszustand nach M1/M2)
- **Mutation:** keine; fatou_fork.gp == fatou.gp
- **Gate:** PASS (Kalibrierlauf, siehe research/gate.py REQUIRED_DIGITS)
- **Benchmark:** bench/results/: baseline-oneshot (warm 355-1091 d/s),
  m1-persistent (warm 42k-136k d/s, Median-Speedup 182x),
  m2-cache (cold 0.05-0.10s dank State-Cache; vorher 1.6-15s)
- **Entscheidung:** Referenzpunkt
- **Learnings:** Speedup-Treppe One-Shot -> persistent -> Cache dokumentiert.
  Pool (M2b): 3.6x auf grossen warmen Batches (4000 Evals), unterhalb ~1000
  Evals frisst der Spawn-Overhead den Gewinn.

---

## exp-001 (2026-07-10) — ctr 8/10 -> 7/10 (Konturterm-Dichte)
- **Hypothese:** lctr (Zeile ~568) skaliert mit 1/(-log(ctr)); ctr=7/10 wird
  durch den ir-Auto-Shrink zu effektiv 0.729 (vs. 0.833 default) ->
  ~42% weniger Konturterme -> Init- und Eval-Pfad deutlich schneller.
  limitp=16 wurde verworfen (kappt Praezision global auf 16 Stellen, reisst
  das Gate sicher). Risiko: groessere Terme/schlechtere Konditionierung ->
  Genauigkeitsverlust; genau dafuer ist das Gate da.
- **Mutation:** fatou_fork.gp Zeile 13: `ctr = 8/10;` -> `ctr = 7/10;`
- **Gate:** PASS -- aber alle Werte IDENTISCH zum Original (Roundtrip auf 6
  signifikante Stellen gleich) = No-Op.
- **Benchmark:** entfaellt (No-Op)
- **Entscheidung:** revert
- **Learnings:** Top-Level `ctr` ist toter Code; sexpinit ueberschreibt es in
  Zeile 425 basisabhaengig (`ctr=8/9` bzw. `ctr=4/5`). Es gibt einen sauberen
  Override-Mechanismus: `myctr`/`myir` (Zeilen 426-427).

---

## exp-002 (2026-07-10) — ctr-Auswahl Zeile 425: 4/5 -> 7/10
- **Hypothese:** weniger Konturterme (lctr ~ 1/-log(ctr)) -> schneller.
- **Mutation:** Zeile 425 else-Zweig `ctr=4/5` -> `ctr=7/10`
- **Gate:** FAIL (28-60 digits statt 74+; nur e|80 fiel bloss 62.7->59.3)
- **Benchmark:** entfaellt
- **Entscheidung:** revert
- **Learnings:** (1) ctr bestimmt das erreichbare Praezisions-Plateau, nicht
  nur Speed. (2) Der e|80-Ausreisser ist fast ctr-unabhaengig -> andere
  Ursache, separat untersuchen. (3) Roundtrip blieb IDENTISCH (2.56296e-87):
  invabel ist Newton-Inversion derselben abel-Funktion -> Roundtrip misst nur
  Newton-Toleranz, NICHT echte Genauigkeit. Agreement vs. Referenz ist der
  einzige echte Genauigkeits-Check.

---

## exp-003 (2026-07-10) — ctr=8/9 fuer alle Basen (Genauigkeits-Versuch)
- **Hypothese:** 8/9 ist der "high-accuracy"-Zweig -> mehr digits (Keep-Regel b).
- **Mutation:** Zeile 425 else-Zweig `ctr=4/5` -> `ctr=8/9`
- **Gate:** FAIL (23-50 digits -- noch schlechter als exp-002)
- **Benchmark:** entfaellt
- **Entscheidung:** revert
- **Learnings:** ir (Z.424) und ctr (Z.425) sind GEMEINSAM abgestimmte Paare
  ((15/16, 8/9) vs. (7/10, 4/5)); einseitige Aenderung zerstoert die
  Kontur-Geometrie (ircircr=ir*ctr*circr). Kontur-Radius-Knoepfe einzeln =
  Sackgasse; Auto-Auswahl ist lokales Optimum. Naechste Familien:
  (A) Init-Staffelung im Fork (Iterationszahlen aus precis), (B) e|80-
  Ausreisser-Fix (gezielte Genauigkeit), (C) Deep-Referenzen dps 500/1000.

---

## exp-004 (2026-07-10) — looplim-Floor auf volle Arbeitspraezision
- **Diagnose vorab:** e|80-Schwaeche ist BREIT (fast alle 32 Cases ~62.7),
  nicht ein Ausreisser-Case. nlim irrelevant (40/60: keine Aenderung);
  looplim=100 hebt auf 79.5. Codelektuere loop() Z.1179ff: `looplim` ist ein
  DIGITS-ZIEL (`while (re<looplim ...)`, re = erreichte Dezimalstellen),
  `nlim` der Iterations-Cap — Semantik invers zur Wrapper-Annahme! Der
  Wrapper bestellt mit looplim=dps-20 nur 60 digits; Basis 2/10 ueberspringen
  das Ziel im letzten Iterationsschritt (Glueck), Basis e stoppt praezise
  beim Ziel. Anomalie am Rande: nlim=50+looplim=100 ergab nur 63.9 (?).
- **Mutation:** fatou_fork.gp loop() nach Z.1190:
  `if ((limitp==0) && (looplim < precis-throwp-2), looplim = precis-throwp-2);`
  (Digits-Ziel = volle Arbeitspraezision, ausser limitp fordert Speed-Cap)
- **Gate:** PASS — e|80: 62.7 -> **79.4 digits (+16.7)**; alle anderen
  Gruppen unveraendert (79.5-80.0 / 199.8); Roundtrips unveraendert.
- **Benchmark:** exp-004 (repeat 3) vs. m3-final, warm-Mediane:
  GESAMT +15.2% (Schwelle 6.8%); e|80 +14.4% UND min digits 62.7->70 (Cap);
  e|200 +88%/2|200 +38% (sub-ms-Timings, mit Vorsicht); 10|80/2|80 neutral.
  Cold rep0 (frischer sexpinit) nicht langsamer als Original.
- **Entscheidung:** **KEEP** (Regel a UND b gleichzeitig erfuellt)

---

## Befund-2026-07-10-deep — Original bei dps 1000 / Basis e KAPUTT
- **Messung (baseline-deep, Original):**
  - e|500: init 46s, warm 48-49k d/s, 478.8 digits (ok, knapp unter Ziel 480)
  - 2|500: init 38s, 479.8 digits (ok)
  - 2|1000: init 227s, 956.3 digits (unter Ziel 980)
  - **e|1000: init 259s, nur 65.9 digits (!!)** — praktisch unbrauchbar
- **Wichtig:** e|500 konvergiert mit denselben nlim=30 Iterationen sauber ->
  bei dps 1000 bricht NICHT die Iterationszahl, sondern etwas anderes
  (Kandidaten: ltht=18-Cap, thsamples-Formel, lctr-Matrixgroesse ~10k,
  stiller theta-init-Fehler). Verbose-Diagnose (quietmode=0) noetig, sobald
  die Maschine frei ist.
- **Frontier-Ziel:** e|1000 auf ~980+ digits bringen = groesster einzelner
  Genauigkeits-Hebel im Projekt.
- **Nachtrag (exp-004-deep):** Fork identisch zum Original auf 500/1000
  (digits UND Init-Zeit) -> looplim bindet dort nicht; exp-004 bleibt keep
  (neutral auf deep). Deep-Warm-Timings sind Einzel-Eval (1 Case/Gruppe,
  ~10-40ms) -> fuer Speed-Aussagen dort groessere Batches noetig.
  Verbose-Diagnose e|1020 + e|520 gestartet (quietmode=0 Konvergenz-Trace).
- **Trace e|520 (Original, nlim=30, looplim=500):** Schleife stoppt am
  nlim=30-Cap bei re=63.9 "decimal digits" (linear ~2.1 digits/Iter,
  ctsamples waechst 20->536, thsamples 7->32). ABER Benchmark mass fuer
  e|500 echte 478.8 digits Agreement -> **re ist das KONTUR-Residuum, nicht
  die finale sexp-Genauigkeit** (renormslog/Downstream verhaelt sich anders).
  Zusammen mit der Anomalie (nlim=50+looplim=100 -> 63.9 digits bei dps 80,
  waehrend nlim=30+looplim=100 -> 79.5) ist die (nlim, looplim)->Genauigkeit
  Landschaft NICHT monoton -> empirisch kartieren. Sweep bei dps 500 laeuft.
- **Trace e|1020 == Trace e|520** (identische Trajektorie, Stopp n=30 bei
  re=63.9): Kontur-Phase ist praezisions-agnostisch. Folgerung: e|500 mit
  478.8 Agreement bei Kontur-Residuum 1e-64 => Agreement dps-vs-dps+20 kann
  KORRELIERTE Fehler verdecken (gleiche Trajektorie beidseits); der e|1000-
  Kollaps auf 66 digits muss DOWNSTREAM der Kontur-Phase entstehen
  (Kandidat: Newton-Cap in betterest/invabel bei hoher Praezision).
  Breakpoint-Sonde dps 600-900 laeuft; looplim bei dps 500 inert bestaetigt
  (480 vs 520 -> identisch 478.8).

---

## Befund-2026-07-10-nlim — Arbeitsmodell: nlim ist die echte Praezisions-Drossel
- **Messungen:** (1) Selbst-Paare e|600..e|900 ALLE exakt 65.9 digits;
  (2) Sweep dps 500: nlim=30 -> \"478.8\" (korreliert), nlim=60 -> \"64.2\"
  (dekorreliert die nlim=30-REFERENZ und zeigt DEREN wahren Fehler);
  (3) Traces: Kontur-Phase stoppt bei n=nlim=30 mit ~536 Termen, ~2.1
  digits/Iter fuer Basis e, praezisions-agnostisch.
- **Modell:** Serienlaenge (nlim-gedeckelt) begrenzt die WAHRE Genauigkeit:
  Basis e: ~2.1 digits/Iter * 30 = ~64-66 digits bei JEDEM dps. Hohe
  Agreements (478.8, 79.4) zwischen gleich-konfigurierten Laeufen sind
  KORRELIERTE Fehler (identische Trajektorien). Basis 2: ~32 digits/Iter *
  30 = ~956 -> erklaert exakt 2|1000=956.3. Komplexe Basen konvergieren
  schnell genug fuer volle Praezision bei dps 80.
- **Konsequenz 1 (Genauigkeit):** nlim muss mit dem Digits-Ziel skalieren
  (Basis-abhaengige Rate; fuer e ~dps/2.1 + Marge). Kandidat exp-005.
- **Konsequenz 2 (KRITISCH, Mensch-Entscheid noetig):** research/reference/
  values.json ist fuer BASIS E jenseits ~digit 66 vermutlich falsch (mit
  nlim=30 erzeugt). Referenz-Regeneration mit konvergiertem Engine waere
  noetig -> per program.md ausserhalb des Loops zu entscheiden; wird dem
  Nutzer vorgelegt sobald das Dekorrelation-Experiment bestaetigt.
- **Laufend:** truth(nlim=30/60) via Vergleich gegen nlim=110@540.
- **Oracle-Versuch paulsen_tetration:** zu schwach — weicht ab digit 10 ab
  (1.6463542333.. vs fatou 1.6463542337..; Kouznetsov-Literaturwert stuetzt
  fatou auf ~14 digits). Kein Urteil ueber digit 66 moeglich.
- **Sweep komplett:** nlim=45/60/100 ALLE exakt 64.2 vs nlim-30-Referenz
  (drei Serienlaengen, konsistent) — Modell weiter gestaerkt. Init-Kosten
  ~quadratisch in nlim (30:52s, 45:125s, 60:260s, 100:897s bei dps 500).
- **Strukturelle Folge:** Der EINGEFRORENE Gate blockiert echte
  e-Genauigkeitsverbesserungen: laengere Serie dekorreliert von der
  nlim-30-Referenz -> e|200 wuerde ~64 statt geforderter 194.8 messen ->
  FAIL. Gate/Referenzen sind auf korrelierten Illusionen kalibriert.
  => Referenz-Regeneration (konvergierter Engine) + Rekalibrierung der
  e-Schwellen noetig = MENSCH-ENTSCHEID (program.md). Entscheidungspaket
  wird vorbereitet; bis dahin nur gate-sichere Speed-Experimente.
- **BESTAETIGT (Dekorrelation gegen nlim=110@540):** truth(nlim=30)=64.2,
  truth(nlim=60)=127.3 digits — Vorhersage 2.1*nlim auf 1% getroffen.
  Rate Basis e = ~2.12 echte digits/Iteration. Init ~quadratisch in nlim
  (110@540: 1259s). Echte e-Tetration: 500 digits ~ nlim 240 (~1.7h),
  1000 digits ~ nlim 475 (mehrere Stunden, einmalig dank State-Cache).
- **Entscheidungspaket laeuft:** research/tools/make_reference_v2.py
  erzeugt values_v2_proposed.json (e|80 + e|200, je mit Dekorrelations-
  Verifikation; e|500 per --dps500 spaeter). Eingefrorene Dateien bleiben
  bis zur Freigabe unangetastet.

---

## exp-005 (2026-07-10) — complextaylor=0 fuer reelle Basen
- **Hypothese:** tht2 (zweite Theta-Serie) ist fuer reelle kc redundant;
  Skip halbiert die Theta-Arbeit pro Init-Iteration -> schnellerer Init.
- **Mutation:** loop() vor initsch(kc):
  `if (imag(kc)==0, complextaylor=0, complextaylor=1);`
- **Gate:** PASS — ALLE Werte identisch zu exp-004-Stand (reell wie
  komplex); tht2 fuer reelle Basen bestaetigt redundant.
- **Benchmark:** Init-Median (frischer Cache, 3 Reps): e|80 3.27 vs 3.28s,
  2|80 1.84 vs 1.84s, e|200 8.55 vs 8.54s — NULL Effekt.
- **Entscheidung:** revert (kein Gewinn; einfacherer Code gewinnt)
- **Learnings:** Theta-Serien (thsamples<=32) sind vernachlaessigbar; die
  Init-Kosten dominiert der staylor/Matrix-Pfad (lctr x lctr Solve pro
  Iteration). Speed-Experimente muessen dort ansetzen (z.B. inkrementelle
  Matrix-Updates, ctsamples-Wachstumsstrategie, weniger Iterationen durch
  besseren Startwert).
- **Postmortem-Nachtrag:** Die Engine setzt complextaylor fuer reelle Basen
  SELBST auf 0 (fatou_fork.gp Z.119/122) — exp-005 war redundant zur
  vorhandenen Logik (verifiziert per Injektions-Probe: complextaylor==0
  nach sexpinit auch ohne Mutation). staylor nutzt den Halbierungspfad
  (samples/2 bei konjugierter Symmetrie) fuer reelle Basen also laengst.
  Vorbereitet: research/tools/apply_reference_v2.py (Ein-Kommando-Migration
  nach Freigabe von DECISION_reference_v2.md Option A).

---

## exp-006 (2026-07-10) — Praezisions-Staffelung frueher Iterationen [GEPARKT]
- **Hypothese:** Iteration n traegt nur ~re digits Signal; fruehe Iterationen
  bei reduzierter realprecision (re+Guard) sparen ~30% Init-Zeit.
- **Profil vorab (dps 200, e):** initsch ~0.13s, KEINE teure Post-Phase;
  Iterationen = 100% der Kosten, spaete dominieren (+0.92s bei 536 samples).
  staylor-Extraktion ist naive DFT O(terms*samples); PARI fft() existiert
  (Radix-2), aber ctsamples!=2^k -> Aufrunden verdoppelt sfunc-Sampling,
  frisst den Gewinn. Deshalb Staffelung zuerst.
- **Mutation:** default(realprecision, re+40) pro Iteration (v1), dann
  raten-adaptiv re+40+3*(re-relast) (v2); Restore vor renormslog.
- **Gate:** FAIL (v1: 65-71 auf schnellen Basen; v2: 68-77). Zwei Ursachen:
  (a) schnelle Basen (10-30+ digits/Iter) brauchen groesseren Guard (echt),
  (b) STRUKTURELL: Staffelung aendert Low-Order-Arithmetik -> Dekorrelation
  von der korrelations-kalibrierten Referenz -> e|200 misst die Wahrheit
  (~64-99) statt der Illusion (199.8) -> FAIL unabhaengig vom Guard.
- **Entscheidung:** revert + GEPARKT bis Referenz-v2-Freigabe (danach ist
  die e|200-Schwelle wahrheitsbasiert ~59 und derselbe Code passiert;
  Guard-Tuning fuer schnelle Basen dann separat verifizieren).
- **Learnings:** Das eingefrorene Gate blockiert inzwischen ZWEI Fronten
  (echte e-Genauigkeit UND Low-Order-beruehrende Speed-Experimente). Die
  Referenz-Entscheidung ist der Flaschenhals des gesamten Loops.
  Nebenbefund: e|80 unter Staffelung 80.0 digits (kein Schaden fuer e).

---

## Befund-2026-07-10-wand — ZWEITE Praezisions-Wand jenseits ~dps 300
- **2|1000-v2 INSUFFICIENT:** main dps=1060 nlim=60 (989s) vs verify
  dps=1080 nlim=70 (2507s) stimmen nur auf **79.9 digits** ueberein
  (noetig 1010). Auch Basis 2 kollabiert bei hohem dps — nlim ist dort
  NICHT die Ursache (60/70 nicht bindend bei Rate ~32/Iter).
- **Muster:** e-Selbst-Paare 600-900 konstant 65.9; 2-Paar bei ~1000: 79.9.
  Dekorrelierte e-Verifikationen bis dps 270 sind dagegen ECHT (99.4/241.8).
  => Es gibt eine zweite, nlim-unabhaengige Genauigkeits-Wand zwischen
  ~dps 300 und ~dps 600 (Kandidaten: interne Taylor-Radien/Sample-Zahlen
  in abelest/invabeltaylor/renormslog, ltht-Cap, subst-Pfad).
- **Diskriminator laeuft:** e|500-v2 (dps 560/580). Zusatz-Probe: Basis 2
  dekorreliert bei dps 500.
- **Konsequenz fuers Entscheidungspaket:** v2-Werte fuer 80/200 sind solide;
  500/1000-Tiers ERST nach Lokalisierung/Fix der Wand.

---

## exp-007 (2026-07-10) — nlim = ceil(looplim) statisch [VERWORFEN -> 007b]
- **Gate:** Basis e SPEKTAKULAER: e|80 79.5 (+15.6), e|200 **193.1 (+129
  echte digits)** — der Genauigkeits-Fix wirkt. ABER alle schnellen Basen
  kollabiert (komplex 30-36, 10|80 36.2, 2|80 41.6, 2|200 41.8):
  Ueber-Iteration jenseits des Plateaus destabilisiert den Kontur-Zustand.
  Sheldonisons nlim=30 war der Schutz der schnellen Basen (erklaert auch
  die alte Anomalie nlim=50+looplim=100 -> 63.9).
- **Entscheidung:** revert -> exp-007b (raten-adaptiv: Cap-Anhebung nur wenn
  Rate bei n=6 < 5 digits/Iter; schnelle Basen behalten Caller-Cap).
- **exp-007b: FAIL, byte-identisch zu 007** — auch schnelle Basen starten
  langsam (< 5 digits/Iter frueh), beschleunigen erst spaeter -> n=6-Regel
  feuerte fuer alle. Diskriminator-Messung: dekorrelierte Wahrheit der
  schnellen Basen bei dps 80 = **80.0 digits (voll, echt)** fuer 2,
  0.8+0.4i, 10 -> deren Referenzen sind REAL; der Kollaps unter 007/007b
  war ECHTE Ueber-Iterations-Degradation. -> exp-007c: inkrementelle
  Verlaengerung nur am Cap-Anschlag solange Gesamt-Rate < 5 und re<looplim.
- **exp-007c: FAIL, wieder byte-identisch.** Trace-Aufloesung: Basis 2 hat
  Kontur-Rate ~1.27/Iter (LANGSAMER als e!) -> Rate-Kriterium feuerte auch
  fuer sie. Schluesselbefund **Faktor-2-Muster**: Basis 2 Wahrheit ~ 2x
  Kontur-re (n=30-Exit: Kontur 41.2, echt 80.0); Basis e Wahrheit ~ 1x
  Kontur-re. Erklaert vermutlich auch die "Wand": 2|500-1000 dekorreliert
  79.9 = 2 x (nlim-30-Kontur ~40) — die Wand ist der nlim-Cap durch die
  Faktor-2-Brille, kein separater Mechanismus! ABER: Verlaengerung fuer
  Basis 2 zerstoert echt (41.6 < 80 gemessen gegen 80-true-Referenz).
  thetamode trennt die Basen nicht (alle 1). Warum e-Kontur=Wahrheit und
  sonst 2x, und warum Ueber-Iteration zerstoert: OFFENE Forschungsfrage.
  -> exp-007d: Verlaengerung NUR fuer e-Familie (abs(kc-1)<1e-6).
- **exp-007d: GATE PASS (voll):** e|80 79.5 (echt, +15.6), e|200 **193.1
  (echt, +128.9)**; alle anderen Basen + Roundtrips exakt unveraendert.
  Benchmark (Keep-Regel b, Speed-Check) folgt auf freier Maschine;
  e|400-Sonde gestoppt (durch Faktor-2-Befund obsolet), e|500-v2 laeuft
  weiter (bleibt korrekt: Referenzen werden mit ORIGINAL erzeugt).
- **exp-007d KEEP (gepaarter A/B-Benchmark unter gleicher Last, repeat 2):**
  GESAMT warm +15.2% (Schwelle 6.8%); e|80 warm +35% und min digits
  62.7->70 (Cap); e|200 min digits 64.2->190 (Cap), warm -4% (Rauschen);
  2|200 +44%, 10|80 +9%, 2|80 +3%. Einmal-Init e|200: 20s->497s unter Last
  (danach State-Cache 1.9s) — ehrlicher Preis echter Praezision.
  ZWEITER Keep des Loops; groesster Genauigkeits-Gewinn.

---

## exp-008 (2026-07-10) — Praezisions-Staffelung mit 2x-re-Guard
- **Hypothese:** exp-006 scheiterte, weil der Guard (re+40) die ECHTE
  Zustands-Praezision der Faktor-2-Basen (2x Kontur-re) abschnitt. Guard
  = 2*re+60 respektiert sie; fruehe Iterationen guenstiger -> zielt auf
  den neuen Schmerzpunkt 497s e|200-Init (nach 007d: n~120 Iterationen).
- **Mutation:** default(realprecision, max(48, min(precis, 2*floor(re)+60)))
  pro Iteration; Restore vor renormslog. (Nutzer: weiterlaufen lassen,
  Experimente fortsetzen — e|500-Verify parallel.)
- **Gate:** FAST PASS — alle dps-80-Gruppen sogar auf 80.0 verbessert,
  e|200 haelt 193.1. Einziger FAIL: 2|200 = 84.7 (gefordert 194.8) —
  das ist ein REVEAL, kein Schaden: 84.7 = 2 x 42 (Faktor-2!).
  Nach 2|200-Rekalibrierung: Re-Gate VOLL PASS.
- **Benchmark (A/B vs 007d-Stand):** Kaltinit 2|200 **-52%** (10->4.8s;
  Basis 2 bleibt ganz im Niedrig-re-Regime), e|80 -5%, 10|80 -3%,
  e|200 nur -2% (die teuren spaeten Iterationen laufen eh voll praezise —
  Erwartung widerlegt). Scheinbare Verluste (warm gesamt -6.9%, 2|200
  \"190->84.7\") sind ARTEFAKTE der bewiesen fake'n 2|200-Referenz
  (Digit-Summe faellt durch Wahrheits-Reveal, nicht durch Tempo).
- **Entscheidung:** **KEEP** (Regel b: Gate-digits ueberall gleich oder
  hoeher — dps-80-Gruppen +0.4-0.5 auf 80.0; reale Timings ueberall
  gleich oder besser). Aggregat-Metrik-Einbruch explizit als
  Referenz-Artefakt dokumentiert, kein Metrik-Gaming.
- **Learnings:** Staffelung lohnt dort, wo der ganze Lauf im Niedrig-re-
  Regime bleibt (Faktor-2-Basen); fuer e-Familie bei hohem dps braucht es
  Kostensenkung der SPAETEN Iterationen (ctsamples-Regime) — naechste
  Front. 12h-Marke erreicht; Loop laeuft weiter.

---

## exp-009 (2026-07-10) — ctsamples-Cap 5*re+120 [REVERT]
- **Profil e|200 (418s):** Rate konstant 2.1/Iter; ctsamples 1988 bei
  re 203; letzte 10 Iterationen je ~10.4s; Stall-Exit verbrennt ~3
  Voll-Iterationen.
- **Gate:** FAIL — e-Gruppen auf ~49 digits: die ~1988 Terme sind ECHT
  noetig (~9.8 Terme/digit fuer Basis e; die 4.3er-Schaetzung galt fuers
  falsche ctr-Regime). stopterms overshootet NICHT. Revert.
- **Learnings:** (1) Kein Fett im Sample-Budget. (2) Bei ~2000 Samples
  kostet 2^k-Padding nur +3% -> FFT-Extraktion (O(n^2)->O(n log n),
  ~170x auf dem Block) wird als exp-011 wieder aktiviert.
  (3) exp-010 zuerst: Stall-Exit straffen (Mini-Gewinne <0.1 digit
  kosten sofort nskip-Credit).

---

## exp-010 (2026-07-10) — Stall-Exit-Straffung
- **Mutation:** `if ((n>3) && ((re-relast) < 0.1), nskip--);` nach
  re-Berechnung.
- **Gate:** PASS — alle Gruppen identisch oder besser (e|80 80.0,
  e|200 haelt 193.1).
- **Benchmark (gezielt):** e|200-Kaltinit 427s -> 400s (-6.3%),
  Wert identisch auf 30 Stellen.
- **Entscheidung:** **KEEP** (Regel b: Genauigkeit identisch, Timing besser).
  VIERTER Keep.
- **Naechster Schritt exp-011:** FFT-Extraktion in staylor (2^k-Padding
  +3% bei ~2000 Samples; O(n^2)->O(n log n) auf dem Extraktions-Block).
  Disziplin: erst Standalone-GP-Prototyp (Koeffizienten-Vergleich alt vs.
  FFT auf 1e-30), dann Fork-Integration, dann Gate.
- **Split-Messung (instrumentierter Lauf, spaete Iterationen e|200):**
  Sampling ~4.4-5.6s, Extraktion+Rest ~4-5s pro Iteration (~50/50).
  FFT-Extraktion projiziert e|200-Init 400s -> ~230-250s. thsamples
  waechst auf 90+ bei dps 200 (nicht 32-gedeckelt wie bei dps 80).
- **exp-011 Prototyp WIP** (scratchpad/fft_proto.gp): komplexer Branch
  GP-Typfehler (t_POL in gtos — zu isolieren), reeller Branch
  Formel-Mismatch bei hohen s (~Magnitude -> Phase/Alias-Fehler).
  Naechster Tick: Minimal-Repro, termweiser Vergleich bei m=4.

---

## Befund-2026-07-10-basis2-decke — Basis 2 hat echte ~80-85-digit-Decke
- 2|200-Referenz war korrelations-fake (199.8); echte Basis-2-Genauigkeit:
  80.0@dps80 (dps-gedeckelt), 84.7@200, 79.9@500, 79.9@1000 — Decke
  ~80-85 digits bei JEDEM dps. Extension hebt sie nicht (exp-007: schadet
  bei 80; 2|1000-Paar nlim 60/70: 79.9). Mechanismus offen (vermutlich
  dieselbe Wurzel wie Faktor-2).
- **Delegierte Entscheidung (Nutzer-Direktive):** REQUIRED_DIGITS["2|200"]
  wahrheitsbasiert auf 79.7 rekalibriert (84.7 - 5). Referenzwert-Datei
  unveraendert (die ersten ~84 digits des 2|200-Werts sind echt und
  tragen den Check). e|500/1000- und 2|500-Referenzen bleiben als
  "jenseits der Decke fake" markiert — Deep-Tier-Neubewertung folgt,
  wenn die Decken-Ursache verstanden ist.

---

## MIGRATION 2026-07-10 — Referenz-v2 Option A (Nutzer-Freigabe im Chat)
- 33 Basis-e-Eintraege in values.json durch dekorrelations-verifizierte
  v2-Werte ersetzt (Provenienz in meta.v2_correction).
- gate.py REQUIRED_DIGITS rekalibriert (wahrheitsbasiert): e|80 58.9
  (Fork misst echt 63.9), e|200 59.2 (echt 64.2). Wieder EINGEFROREN.
- Verifikation: Gate PASS, Slow-Gate-Tests (pass+sabotage) gruen,
  Schnell-Suite 48 gruen.
- Genauigkeits-Front ist FREI: Verbesserungen der echten e-Genauigkeit
  (z.B. nlim-Skalierung im Fork) sind jetzt gate-messbar (Keep-Regel b).
  exp-006 (Staffelung) entparkt. Deep-Tiers (500/1000) warten auf
  Wand-Lokalisierung.

---

## Befund-2026-07-10-rate — Konvergenzrate faellt mit n (Modell-Korrektur)
- v2 e|80: verifiziert 99.4 digits (32 Werte) OK.
- v2 e|200 mit nlim=120: nur 193.1 verifiziert (Vorhersage 254) ->
  marginale Rate im Bereich 60->120 nur ~1.1 digits/Iter (statt 2.1).
  Lineare Extrapolation war optimistisch; hohe Ziele (500/1000 digits)
  brauchen ueberproportional mehr Iterationen. Retry e|200 mit 170/250
  laeuft; Rate-Kurve wird aus den Verify-Ergebnissen mitgeschaetzt.
- **KORREKTUR (Retry-Ergebnis):** nlim=170/250 lieferte IDENTISCHE Zeiten
  und identische 193.1 wie 120/180 -> nlim war nicht bindend! Die Schleife
  exitet am looplim=precis-throwp (Kontur-re ~220 bei dps 220). WAHRE
  Genauigkeit = Kontur-re MINUS skalenabhaengiger Gap: ~0 bei re 64,
  ~0.6 bei re 100, ~27 bei re 220. Fuer X echte digits: ARBEITSPRAEZISION
  ~ X + throwp + gap(X) + Marge; nlim nur als Cap >~ re/2 halten.
  Kosten O(precis^2)-ish ueber die Iterationszahl. exp-004-Richtung
  (looplim-Floor) damit voll bestaetigt. Retry 2: main dps 250/verify 270.
- **Learnings:** Wrapper-Bestellung looplim=max(35, dps-20) war die Ursache
  der e|80-Schwaeche. Roundtrip-Identitaet ueber alle Experimente bestaetigt
  erneut: nur Agreement-vs-Referenz misst echte Genauigkeit.
