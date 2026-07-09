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
