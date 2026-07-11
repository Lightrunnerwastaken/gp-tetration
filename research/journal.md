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
- **exp-011/011b (FFT-Extraktion):** Prototyp exakt (2.4e-71). Vollintegration
  (011) zerstoerte Faktor-2-Basen (30-43 digits + log(0)-Crash — der
  2^k-Grid-Wechsel bricht deren Sample/Terms-Koevolution); 011b beschraenkt
  FFT+2^k auf die e-Familie (efam-Global aus loop()). Gate: VOLL PASS auf
  Bestwerten (80.0-Sweep, 84.7, 193.1). A/B gezielt: e|200-Kaltinit
  495s -> 392s (**-21%**), Werte identisch. **KEEP (#5).** Verbleibende
  Luecke zur 45%-Projektion = 2^k-Padding-Kosten im Sampling; Sampling
  (sfunc-Auswertungen) ist jetzt der dominante Block -> exp-012-Kandidat:
  inkrementelle Sample-Updates via Theta-Delta.
- **exp-013-Probe (Basen-Landkarte):** dekorrelierte Wahrheit @dps80:
  1.5/5/7.4 = 80.0 (voll), **Basis 3 = 62.4** — kc=log(log(3))+1~1.094
  liegt NAHE der e-Familie: Langsam-Konvergierer ausserhalb des engen
  007d-Fensters (1e-6). Faktor-2/Decke betrifft weiterhin nur die
  schnellen Basen und erst jenseits ~dps 85.
- **exp-014 (Fenster-Weitung |kc-1|<0.12):** Basis 3: 62.4 -> **80.0**
  echte digits (dekorreliert verifiziert); Fenster schliesst 2 (0.63)
  und 10 (1.83) weiter aus. Gate: VOLL PASS (alle Gruppen Bestwerte).
- **Entscheidung: KEEP (#6)** (Regel b: Basis 3 +17.6 echte digits, keine
  Verluste; Standard-Gruppen unberuehrt).
- **kc-Landkarte komplett (Fork exp-014, dekorreliert @dps80):**
  1.5/1.7/2.5/4/5/7.4 alle 80.0 echte digits; 2.5 (kc=0.913, IM Fenster)
  unbeschadet extended; 4 (1.327) und 1.7 (0.366) brauchen keine Extension.
  Fenster |kc-1|<0.12 sitzt richtig: Langsam-Konvergierer e (kc=1) und
  3 (1.094) drin, gesunde Basen draussen/unbeschadet.
- **exp-012 (sfunc-Walk-Cache) implementiert:** Der fs/finv-Walk in sfunc
  ist ct/theta-UNABHAENGIG (nur Basis-Map) -> Endpunkte (z, n) pro
  Sample-Index cachebar, solange das Grid unveraendert ist (e-Familie:
  2^k-Plateaus!). Globals swon/swidx/swkey/swz/swn/swvalid; Cache nur im
  staylor-Sampling der e-Familie aktiv. Smoke: e und 2 wertidentisch.
  Gate laeuft. Erwartung: ~Haelfte des Sampling-Blocks in Plateau-
  Iterationen gespart.
- **Keep-Treppe exp-012/015/016 (sfunc-Redundanzen, alle Gate-PASS auf
  Bestwerten, Werte identisch):** 012 Walk-Cache -6.4% (359->336s),
  015 lazy y2 -13.2% (340->295s), 016 Base-Est-Cache -4.3% (299->286s).
  **Kumulativ e|200-Kaltinit seit 007d: 497 -> 286s (-42%) bei +129
  echten digits.** sfunc-Kette abgeerntet (verbleibend nur der
  ct-abhaengige Inside-abelest — nicht cachebar). Naechste grosse Fronten:
  Faktor-2-Wurzel (Genauigkeit jenseits ~85 fuer schnelle Basen),
  Deep-Tier-Referenzen (jetzt ~40% billiger).
- **Hypothesen-Kette Faktor-2/Decke (naechster Forschungs-Tick):**
  (a) re misst die UPDATE-Groesse der Iteration, nicht den Fehler; bei
  quasi-linearer Konvergenz (e-Familie, Kontraktion nahe 1) ist
  Fehler ~ Update (Faktor 1x), bei staerkerer Kontraktion waere
  Fehler ~ Update^2 (2x digits) — passt zu Faktor-2, ABER erklaert die
  starre ~80-85-Decke von Basis 2 nicht (mehr Iterationen muessten sie
  quadratisch heben; beobachtet: 79.9 konstant, Extension schadet sogar).
  (b) Deshalb Verdacht: Die Decke kommt aus der EINMALIGEN
  initsch-Matrixloesung (m2 x = s2, lctr x lctr) — deren Loesung koennte
  bei ~80-85 digits Genauigkeit liegen (Vandermonde-artige
  Konditionierung!) und alles Downstream deckeln; e-Familie evtl. besser
  konditioniert (anderes ircircr-Regime). TEST: Residuum ||m2 x - s2||
  nach dem matsolve direkt ausgeben (Basen 2 vs e, dps 200) — wenn
  ~1e-84 fuer Basis 2 -> Decke lokalisiert; Fix-Kandidat: iterative
  Nachverfeinerung der Loesung (residual correction, 1-2 Schritte).
  KORREKTUR nach Code-Check: der grosse lctr-Solve (Z.690) liegt auf dem
  matrixradius-SONDERPFAD (Prints erschienen in keinem Standard-Trace);
  aktive matsolves (Z.489/540) sind 2x2-trivial. Hypothese (b) faellt.
  Decken-Kandidaten offen: Theta-Serien-Genauigkeit (thsamples-Formel),
  ircircr-Geometrie, superf-Seed-Genauigkeit. Naechster Test: thsamples
  bei Basis 2 kuenstlich erhoehen (tht_re_mult-Faktor) und Decke messen.
- **DECKE GEBROCHEN (Geometrie!):** Ausschluss-Serie: thsamples x2.5 -> 83.2
  (kein Effekt), superfk=16 -> 83.2 (kein Effekt). ABER konservatives
  Kontur-Paar via initsch(kc, 8/9, 15/16): Basis 2 dekorreliert @200 =
  **113.9 digits** (Decke war 84.7, +29!). Die Basis-2-Decke ist
  GEOMETRISCH (Konturradius vs. Singularitaeten des theta-Mappings).
  Eskalations-Test (19/20, 24/25) mit nlim=120 laeuft. Trade-off beachten:
  dichteres ctr = mehr Terme/digit (bis 19.5 statt 4.3) = teurere Serien.
  Fix-Kandidat exp-017: geometrie-adaptive Paar-Wahl nach Ziel-dps.
  Eskalation (19/20, 24/25): FEHLSCHLAG — 0s-Inits + nur 66.7 digits
  (lctr ~45 Terme/digit -> struktureller Bruch mit stillem Fallback).
  exp-017-Design also: Standard-Paar fuer dps<=~90, (15/16, 8/9) fuer
  hoehere Ziele bei schnellen Basen; Kosten-Check noetig (19.5 statt
  4.3 Terme/digit -> Serien ~4.5x laenger, warm-Eval mit betroffen).
- **exp-017 (Eskalation precis>120): Gate-FAIL 2|200=38.9** — Gate ruft
  mit nlim=30; dichte Geometrie konvergiert langsamer -> Cap wuergt ab.
- **exp-017b (Eskalation + Extension gekoppelt): Gate-FAIL 2|200=41.8**
  (~exp-007-Zahl!) trotz teurem Lauf. Verdacht: exp-010-Stall-Straffung
  ((re-relast)<0.1 -> nskip--) toetet die dichte Geometrie in ihrer
  langsamen/oszillierenden Anfangsphase (Probe mit hartem nlim=60 lief
  durch -> 113.9). Verbose-Trace des Gate-Setups laeuft zur Verifikation;
  Fix-Kandidat 017c: Stall-Straffung bei geoesc lockern (Schwelle 0.01)
  oder nskip-Startkredit erhoehen.
- **AUFLOESUNG (Trace): Stall-These FALSCH, 113.9-Probe war ILLUSION.**
  Trace: dichte Init laeuft sauber bis Kontur-re=199 (n=168, Rate ~1.2)
  — und das Gate misst trotzdem 41.8. Der Gate-Vergleich (dichte vs.
  alte Geometrie) ist bereits KREUZ-dekorreliert -> dichte Geometrie
  liefert real ~42. Die 113.9-Probe verglich zwei Laeufe DERSELBEN
  Geometrie (nur dps verschieden) -> korrelierte Eval-Fehler.
  VERSCHAERFTE LEKTION: Dekorrelation braucht METHODEN-Diversitaet;
  dps-Diversitaet allein reicht nicht.
- **exp-017/b: REVERT.** Wichtigster Neben-Ertrag: Der ~84-Deckel der
  schnellen Basen sitzt im EVAL-PFAD (sexp->invabel->betterest/abelest,
  ircircr-Radien), NICHT in der Kontur-Serie (Kontur-re 199 hilft
  nichts). Naechste Decken-Front: Eval-Pfad-Analyse (Newton-Toleranzen,
  Radien-Entscheidungen, betterest-Konvergenz bei hohem dps).
- **Eval-Pfad-Analyse (Faktor-2-Wurzel, Stand):** Messungen Basis 2@200:
  invabel-Evalpunkt bei 0.405*circr (INNERHALB Sampling 0.56) — Serien-
  Extrapolation als Deckel widerlegt; Serie nur 304 Terme (nlim-30-Exit).
  KOMPENSATIONS-THESE: Faktor-2 = Fehlerkompensation 1. Ordnung — sexp =
  invabel(z - rslog), rslog wird aus DERSELBEN fehlerhaften Abel-Funktion
  kalibriert (renormslog(ct)) -> korrelierte Fehler heben sich linear,
  Rest O(err^2) = doppelte digits. Direkter rslog-Stoertest zeigt
  Sensitivitaet ~1 (widerlegt nichts: Stoerung war UNkorreliert).
  RICHTIGER TEST (naechster Tick): ct-Serie stoeren (z.B. +1e-50*x^2),
  rslog NEU via renormslog(ct) kalibrieren, dann delta_sexp messen —
  wenn <<1e-50: Kompensation bestaetigt; Konsequenz waere, dass ECHTE
  200 digits fuer schnelle Basen eine sauber konvergierte Kontur brauchen
  (Extension ohne Degradation — deren Ursache dann die letzte Frage).
  ERGEBNIS: Daempfung nur ~45x (2.2e-52) — KEINE quadratische
  Ausloeschung; einfache Kompensations-These NICHT bestaetigt.
  Verfeinerung: echte Kontur-Fehler sind strukturierte Fourier-Reste,
  kein glattes x^2 — sauberer Test braeuchte einen echten Fehler-Vektor
  (Differenz zweier ct-Staende milder nlim-Differenz, z.B. ct_28 vs
  ct_30, mit jeweils passendem rslog). Faktor-2-Wurzel bleibt OFFEN;
  alle Messwerkzeuge und Ausschluss-Ergebnisse dokumentiert.
- **DURCHBRUCH (Fehler-Vektor-Test): sexp-Fehler = Kontur-Fehler, FAKTOR 1.**
  delta_sexp(nlim28 vs 30, Basis 2@200) = 7.9e-40 == Kontur-Niveau des
  schlechteren Standes (re_28=38.7). KASKADEN-KONSEQUENZEN:
  (1) 2|200 ist real nur ~41 digits genau — die 84.7-Messung teilte den
  SYSTEMATISCHEN Trunkierungsfehler beider nlim-30-Seiten und mass nur
  Rundungs-Divergenz (Illusion Nr. 3a). Faktor-2-Muster war Artefakt.
  (2) Die "Ueber-Iterations-Degradation" (exp-007: 41.6) war Illusion 3b:
  Extension lieferte Kontur~199-Werte, die die ~41-true-REFERENZ
  entlarvten. Extension schadet NIE — exp-007-statisch (alle Basen) war
  die RICHTIGE Mutation, abgelehnt am fehlerhaften Massstab!
  (3) dps-80-Werte bleiben gueltig (looplim-Exit: Kontur 79.5 = echte 80).
  (4) Dekorrelations-Regel final: NUR der Fehler-Vektor-Test (Differenz
  echter Iterations-Staende) misst Wahrheit; dps-Paare UND Arithmetik-
  Divergenz teilen die Systematik.
- **PLAN (delegierte Befugnis):** (a) 2|200-Referenz neu: Original mit
  nlim hoch + looplim=0 @dps250, Konvergenz via Stand-Differenz-Test;
  (b) values.json ersetzen + 2|200-Schwelle wahrheitsbasiert anheben;
  (c) exp-007-statisch (Extension fuer alle Basen) erneut gaten -> der
  Weg zu echten 200+ digits fuer ALLE Basen ist frei.
- **exp-018-Gate: 2|200 = 194.3 PASS gegen die bewiesene Referenz** —
  Extension-fuer-alle experimentell validiert. GLEICHZEITIG Illusion 3c:
  dps-80-Gruppen der schnellen Basen FAILen mit Reveal-Zahlen (30-42),
  weil auch bei dps 80 der nlim=30-Cap VOR looplim exitet (Basis 2
  braucht n~63 fuer Kontur 79.5) -> ALLE dps-80-Referenzen der schnellen
  Basen sind real nur ~30-42 digits; die 80.0-Dekorrelations-Pruefungen
  teilten die nlim-30-Systematik (dps-Paare beweisen nichts — final
  bestaetigt). SANIERUNG v4: alle schnellen-Basen-Referenzen konvergiert
  (looplim=0, nlim non-binding) + Fehler-Vektor-verifiziert neu erzeugen;
  bei ref_dps 100 billig. Danach Schwellen final wahrheitsbasiert.
- **FINAL-GATE VOLL PASS (alles echt): exp-018 KEEP (#10).**
  Alle dps-80-Gruppen 80.0 echte digits (gegen 99.4-100.0-bewiesene
  v4-Referenzen); 2|200 = 194.3 echt (heute frueh real ~41!); e|200 193.1.
  Das Mess-System (Referenzen + Gate) ist erstmals durchgehend
  wahrheitsbasiert; Extension gilt fuer alle Basen. Ehrlicher Preis:
  Gate-Lauf jetzt ~25-40 min (echte Konvergenz aller Gruppen), Inits
  einmalig teurer — State-Cache amortisiert.
- **e|500-v2 GESTORBEN am 4h-init_timeout** (dps 560, nlim 340 braucht
  >4h/Lauf auf dieser Maschine); der blinde WorkerDied-Retry hat den
  Timeout VERDOPPELT (~8h verbrannt). Fixes: (1) Wrapper retryt
  Init-Phase-Tode nicht mehr (nur Eval-Phase), (2) Deep-Tier-Referenzen
  ZURUECKGESTELLT bis FFT-Extraktion die Init-Kosten ~halbiert; dann
  Neuplanung (ggf. Referenz-Policy-Amendment: fork-erzeugte Referenz mit
  dekorrelierter Original-Kreuzverifikation).

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

---

## Referenz-v5 (2026-07-11) — Erste bewiesene e|500-Referenz (Marathon-Ergebnis)
- Original-Engine-Marathon (sequentiell dps 520 + dps 533, ~16h Wandzeit)
  abgeschlossen: kontur_e533 = 514.8.
- **Fehler-Vektor (520 vs 533): 1.748e-497** -> der dps-533-Wert
  sexp_e(0.5) hat >= 497 bewiesene digits. Erster bewiesener Deep-Tier-Wert.
- Alter v1-Eintrag sexp|e|0.5|500 divergierte ab Digit 63 (gleiche
  ~64-digit-Systematik wie die alte e|200-Referenz) -> ersetzt,
  Provenienz in meta.v5_correction (delegierte Entscheidung).
- Noch offen im 4-Wege-Rennen: Original 2|500 (dps 520+533), Fork-Racer
  e|520 + 2|520 (Engine-Diversitaet als Zusatz-Beweis; Original laeuft
  weiter, nicht abgebrochen — Nutzer-Direktive).

---

## Kalibrierung-2026-07-11 — Fork Digits-vs-dps-Kurve gegen bewiesene e|500-Referenz
- Probe dps 300 (Fork, sexpinit(e,400,4,0), unter Last von 4 parallelen
  gp-Jobs): **294 echte digits** (Vergleich gegen v5-Referenz >=497 bewiesen),
  Laufzeit ~15-20 min. looplim-Ziel (~280) um 14 digits ueberschossen ->
  Gap-Modell im Deep-Bereich guenstiger als bei dps 220 (dort 193 echt).
  Faustregel neu: echte digits ~ dps - 6 (statt dps - 27) bei dps 300.
- Probe dps 400 laeuft (Erwartung ~390+ echte digits).

---

## Engine-Diversitaet-2026-07-11 — Fork bei 500 digits korrekt (staerkste Beweisform)
- Fork-Racer e|dps520 (alle 10 Keeps) fertig nach ~5.4h: kontur_re 496.39.
- **Cross-Check vs Original-Referenz (dps 533, >=497 bewiesen): 495 digits
  Uebereinstimmung.** Zwei verschiedene Engines, verschiedene Tiefen ->
  e|500-Referenz engine-divers bestaetigt; Keep-Stack korrumpiert Deep-Tier
  NICHT (FFT-Extraktion, Staffelung, Caches alle sauber bei dps 520).
- Faktor-1-Modell haelt im Deep-Bereich: Kontur 496.3 ~ gemessene 495.
- Speed-Indiz: Fork ~5.4h vs Original ~8h/Lauf (beide Mischlast) -> ~1.4x.
- In meta.v5_correction.engine_diversity_check dokumentiert.

- Probe dps 400: **382 echte digits** (~83 min unter Last). Kurve komplett:
  220->193, 300->294, 400->382, 520->495. Modell: echte digits ~ dps - 24
  + Ueberschuss 0..+20 (Iterations-Quantisierung; Exit bei re>=looplim,
  looplim ~ dps-24). Fuer X garantierte digits: dps ~ X+30 bestellen.
- Zeitskalierung Deep-Tier steil: 18min/83min/5.4h fuer dps 300/400/520
  (Mischlast) -> Speed-Hebel ist Zeit pro Iteration (Serie/sfunc), nicht
  Praezisions-Bestellung. Naechster Angriffspunkt nach 2|500-Abschluss.

---

## Profiling-2026-07-11 — Deep-Tier-Zeitverteilung (research/tools/fatou_profile.gp)
- Instrumentierte Fork-Kopie, dps 300 (1606s Last) + dps 220 (Feinprofil).
- **staylor = 92% der Gesamtzeit; davon ~98% die sfunc-Sampling-Schleife.**
  FFT-Extraktion (exp-011b) ist erledigt (~0.2s/Iter vs 12s Sampling).
- stopterms waechst linear: ~9 Terme/Digit (re 105 -> st 933, re 205 ->
  st 1945). Power-of-2-Rundung (efam-FFT-Grid) verschwendet bis ~2x:
  re 217 braucht ~2170 Terme -> Grid 4096 (2048 sfunc-Calls statt ~1086).
  Ratchet bestaetigt: length(ct) zwingt Folge-Iterationen aufs grosse Grid;
  letzte Iteration dps300 lief auf Grid 8192.
- Konvergenz strikt linear ~2.05 digits/Iter (144 Iter fuer 294 digits).
- Instrumentierung validiert: Profil-Lauf liefert identische 294 echte
  digits wie saubere Probe.
- Abgeleitete Experimente: exp-019 Aitken-Delta^2 auf ct-Koeffizienten
  (Iterationen ueberspringen), exp-020 exakte Grids (Bluestein) gegen
  Power-of-2-Verlust. Fernziel-Hebel: Newton-Kantorovich (digits verdoppeln
  statt +2.1/Iter) und Residuen-Auswertung (Praezisions-Leiter).

---

## exp-019 (2026-07-11) — Aitken-Delta^2 auf ct-Koeffizienten: REVERT
- Mutation: alle 6 Iterationen Aitken-Extrapolation ueber die letzten drei
  ct-Iterierten (Fenster re in [40, looplim-10]), rr danach neu berechnet.
- A/B dps 150 (paarweise, gleiche Last): aitk0 81.7s/72 Iter vs aitk1
  77.9s/71 Iter, beide ~147 echte digits -> 4.6% < Rauschschwelle, REVERT.
- **Mechanismus-Erkenntnis (wertvoll):** Die Kontur-Iteration ist KEINE
  Kontraktion mit dominantem Mode in festem Raum, sondern eine
  informationslimitierte Serien-Erweiterung: jede Iteration erzeugt ~2
  digits durch NEUES Sampling (st waechst ~9 Terme/digit mit). Extrapolation
  kann ungesampelte Information nicht erzeugen -> Iterationszahl ist durch
  den Informationsfluss gepinnt. Beschleunigung => Sample-Pass verbilligen
  (Grid exakt statt Power-of-2 = exp-020; sfunc-Kosten) oder Newton-artige
  Schritte mit mehr Informationsgewinn pro Pass (Forschung).

---

## exp-020 (2026-07-11) — 256er-Grids + Bluestein-Extraktion: KEEP (#11)
- Mutation (efam-only): staylor-Grid > 256 wird auf 256er-Vielfache
  quantisiert statt Power-of-2 (Verschwendung <=12% statt bis 2x im
  dominanten sfunc-Sampling); Extraktion routet bei Nicht-Power-of-2
  auf bluedft() (Bluestein-Chirp + PARI-Poly-Mult, exakt-N-DFT).
  Prototyp research/tools/bluestein_proto.gp: relerr ~1e-73 @dps60,
  53x schneller als Rotations-DFT bei N=2100.
- A/B gepaart (identische Last): dps 220: 461.2s -> 353.9s = **1.30x**
  (Grid 8192 -> 2560), beide 191 echte digits. dps 150: 81.7 -> 68.3s =
  **+16.4%**. Gewinn waechst mit Tiefe (Grid-Waste waechst).
- Gate: FULL PASS (alle 80.0/194.3/193.1, Roundtrips sauber). Suite 45
  passed, 3 skipped. exp-012-Walk-Cache bleibt innerhalb der 256er-Stufen
  gueltig (Grid-Identitaet [samples,w,r] unveraendert pro Stufe).
- Erwartung Deep-Tier: ~1.4x bei dps 300+; naechster Hebel = sfunc-Kosten
  pro Call bzw. Newton-artige Schritte (mehr Info pro Sample-Pass).
