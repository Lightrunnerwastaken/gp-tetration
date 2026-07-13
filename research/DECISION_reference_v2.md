# [ENTSCHIEDEN 2026-07-10: Option A — vom Nutzer im Chat freigegeben und umgesetzt]
# Entscheidung nötig: Basis-e-Referenzwerte sind jenseits ~digit 64 falsch

Stand: 2026-07-10, Autoresearch-Loop (Details + Messreihen: `research/journal.md`)

## Was passiert ist (Beweiskette, jeweils unabhängig gemessen)

1. **Selbst-Paare** (Engine bei dps X vs. dps X+20, Basis e): 600–900 → alle
   exakt 65.9 digits Übereinstimmung. Bei dps 500 dagegen scheinbar 478.8.
2. **Dekorrelation**: Läufe mit nlim=45/60/100 gegen die nlim=30-Referenz →
   alle exakt 64.2. Gegeneinander (nlim 30 vs. 110): truth(30)=64.2,
   truth(60)=127.3 — die Vorhersage des Modells (2.1·nlim) auf 1 % getroffen.
3. **Ursache**: Der Wrapper übergibt überall nlim=30 (Iterations-Cap).
   Basis e konvergiert mit ~2.1 echten digits/Iteration → jede mit nlim=30
   erzeugte Basis-e-Zahl hat nur ~64 echte digits — auch die eingefrorenen
   Referenzwerte (nominell 100–1020 digits). Hohe gemessene Agreements waren
   korrelierte Fehler (identische Iterations-Trajektorien beidseits).
   Basis 2 ist milder betroffen (~32 digits/Iter → 2|1000 real ~956),
   komplexe Basen konvergieren schnell genug (nicht betroffen bei dps 80).

## Warum das nicht im Loop fixbar ist

`research/gate.py` + `research/reference/values.json` sind per program.md
eingefroren (Anti-Gaming). Aber: Eine ECHTE Genauigkeitsverbesserung des Forks
dekorreliert von der falschen Referenz und würde z. B. bei e|200 nur ~64
messbare digits zeigen → Gate-FAIL trotz besserer Wahrheit. **Die falsche
Referenz blockiert die Genauigkeits-Front strukturell.**

## Vorschlag (Dateien liegen bereit, nichts Eingefrorenes wurde angefasst)

`research/reference/values_v2_proposed.json` enthält korrigierte Basis-e-Werte
mit konvergiertem Engine, jeder Wert per Dekorrelation verifiziert:

- 32× e|80-Cases: verifiziert auf **99.4** digits (nötig: 90) ✓
- 1× e|200-Case: verifiziert auf **241.8** digits (nötig: 210) ✓
- e|500 rechnet gerade (korrigiertes Verfahren, main dps 560)
- e|1000: nach gleichem Muster, Kostenpunkt mehrere Stunden (ein Lauf)

**Bei Freigabe würde ich:**
1. Die Basis-e-Einträge in `values.json` durch die v2-Werte ersetzen
   (Provenienz-Meta dokumentiert die Verifikation).
2. `gate.py` REQUIRED_DIGITS für e-Gruppen neu kalibrieren (der jetzige
   Fork misst gegen WAHRE Referenzen: e|80 ~79 → Schwelle ~74;
   e|200 ~193 → Schwelle ~188) — einmalig, dokumentiert, dann wieder
   eingefroren.
3. Benchmark-Historie ab dann mit v2-Referenzen fahren (Label-Konvention
   markiert den Schnitt).

## Optionen

- **A (empfohlen): Freigeben** — Referenzen v2 übernehmen + Gate-Schwellen
  rekalibrieren. Genauigkeits-Front wird frei; „mehr echte Stellen" wird
  erstmals messbar.
- **B: Nur additiv** — v2 als Zusatz-Gate führen, altes Gate unverändert
  (Loop müsste beide bestehen; e-Genauigkeitsverbesserungen bleiben am
  alten e|200-Kriterium blockiert — nicht empfohlen).
- **C: Ablehnen** — Loop arbeitet nur an der Speed-Front weiter.

Antwort im Chat genügt („Option A" o. ä.); bis dahin läuft der Loop
gate-sicher weiter (Speed-Experimente + v2-Vervollständigung).
