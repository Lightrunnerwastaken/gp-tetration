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
