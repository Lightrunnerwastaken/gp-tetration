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
- **Learnings:** Wrapper-Bestellung looplim=max(35, dps-20) war die Ursache
  der e|80-Schwaeche. Roundtrip-Identitaet ueber alle Experimente bestaetigt
  erneut: nur Agreement-vs-Referenz misst echte Genauigkeit.
