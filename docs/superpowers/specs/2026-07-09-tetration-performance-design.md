# Design: Tetration-Performance maximieren (Stufe 1 Engineering + Stufe 2 Autoresearch-Loop)

Datum: 2026-07-09
Status: Entwurf zur Review

## Ziel

Mehr verifizierte Dezimalstellen pro Sekunde aus dem fatou.gp-Rechenkern holen — sowohl
Durchsatz (viele Auswertungen bei ~80–150 dps, mixed_phase-Workload) als auch Präzision
(Einzelwerte bei 200–1000+ dps). Messgröße für alles: **verifizierte korrekte Stellen
pro Sekunde** auf einem festen Workload-Mix.

## Getroffene Entscheidungen

- **Beides optimieren** (Durchsatz und Präzision), erst strukturelle Bremsen, dann Kern.
- **fatou.gp darf geforkt und mutiert werden**; das Original bleibt eingefroren als
  Korrektheits-Referenz.
- **Benchmark reell-fokussiert** (Basen e, 2, 10); komplexe Basen (1+I, 0.8+0.4i, 2+I)
  bleiben als Korrektheits-Anker im Gate.
- **Gestufter Hybrid**: Stufe 1 deterministisches Engineering (kein Harness), Stufe 2
  Autoresearch-Loop im Karpathy-Muster (karpathy/autoresearch als Vorbild, nicht als
  Codebasis — dessen Code ist GPU/GPT-spezifisch).
- **Loop-Fahrer ist Claude** direkt im Repo (optional autonom via /loop). Die eigene
  autoresearch-Harness (`C:\Users\janis\autoresearch`) kann später als Fahrer andocken,
  ist aber kein Blocker.

## Stufe 1 — Engineering (sichere 10–100×)

### 1.1 64-bit PARI

`find_default_gp_exe()` bevorzugt heute Pari32 (`gp_backend.py`). Reihenfolge drehen
(Pari64 zuerst), Pari64 installieren. 32-bit-Limbs kosten bei hoher Präzision grob
Faktor 2–4. Kein API-Bruch.

### 1.2 Persistente GP-Session (`FatouGPWorker`)

Heute startet jeder Batch `gp.exe` neu und läuft `sexpinit` komplett neu; bei dps 80
dominiert die Initialisierung die Gesamtzeit. Neu:

- Ein langlebiger `gp.exe`-Prozess pro (Basis, dps) via `subprocess.Popen` mit
  stdin/stdout-Pipes.
- Protokoll wie heute markerbasiert (`__RES__n`, real/imag als getrennte Zeilen),
  nur ohne Prozess-Neustart; `sexpinit` läuft genau einmal beim Spawn.
- `FatouGP`-API bleibt unverändert; `FatouGPSession` nutzt intern den persistenten
  Worker. `mixed_phase` und CLI profitieren ohne Umbau.

Fehlerbehandlung: Timeout pro Auswertung; bei Hänger oder Parse-Fehler Worker killen,
neu spawnen, Fehler mit GP-Rohausgabe nach oben. Prozess-Lebensdauer an das
Python-Objekt gebunden (Context-Manager + `__del__`-Fallback).

### 1.3 sexpinit-State-Cache

- Nach erfolgreichem `sexpinit` die State-Globals (Theta-/Taylor-Koeffizienten) per
  PARI `writebin` nach `~/.cache/fatou_backend/` schreiben.
- Cache-Schlüssel: Hash aus (Basis-Ausdruck, dps, nlim, nskip, looplim,
  fatou.gp-Datei-Hash).
- Beim Start: `readbin` statt Neuberechnung, danach Selbst-Validierung gegen einen
  mitgecachten Ankerwert (`sexp(0.5)`). Validierung fehlgeschlagen → Cache verwerfen,
  voll initialisieren.
- Offener Punkt (empirisch in der Implementierung zu klären): welche GP-Globals bilden
  den vollständigen sexpinit-State. Vorgehen: Init → Dump aller Top-Level-Variablen →
  Restore in frischem Prozess → Ergebnisvergleich.

### 1.4 Worker-Pool

`FatouGPPool(base, n_workers)` verteilt Batches round-robin auf N persistente Worker.
Zielnutzer: `mixed_phase` mit hohen Sample-Zahlen.

### 1.5 Benchmark (wird von Stufe 2 unverändert geerbt)

`bench/benchmark.py` mit festem Workload-Mix:

- Durchsatz: sexp/slog-Batches, dps 80, Basen e/2/10.
- Präzision: Einzelwerte, dps 200/500/1000, Basen e/2.
- Komplexe Anker: nur Korrektheits-Check, fließen nicht in die Speed-Metrik ein.

Metrik: verifizierte korrekte Stellen pro Sekunde, verifiziert gegen einmalig erzeugte
Referenzwerte der eingefrorenen Original-fatou.gp bei höherem dps (dps+20).
Kaltstart (inkl. Init) und Warmpfad (persistente Session) werden getrennt ausgewiesen.

### Tests Stufe 1

- Bestehende 12 Tests bleiben grün (API unverändert).
- Neu: Worker-Lebenszyklus (Spawn/Timeout/Restart), Cache-Hit/-Miss/-Korruption,
  Pool-Äquivalenz (Pool-Ergebnisse == Einzel-Worker-Ergebnisse).

## Stufe 2 — Autoresearch-Loop (Karpathy-Muster)

### Struktur

```
research/
  program.md      # Ziel, Regeln, Mutations-Leiter — Betriebsanleitung des Loops
  gate.py         # UNVERÄNDERLICHES Korrektheits-Gate
  reference/      # eingefrorene Referenzwerte (Original-fatou.gp, dps+20, eingecheckt)
  journal.md      # Experiment-Log: Hypothese → Mutation → Gate → Metrik → keep/revert
bench/benchmark.py  # aus Stufe 1
src/fatou_backend/vendor/
  fatou.gp        # ORIGINAL, eingefroren (Referenz)
  fatou_fork.gp   # Mutationsgegenstand, startet als Kopie
```

### Experiment-Schleife

Ein Experiment = genau eine Mutation →

1. `gate.py`: Übereinstimmung mit Referenzwerten bis zur geforderten Stellenzahl,
   Roundtrip-Residuen unter Schwelle, komplexe Anker korrekt.
2. Nur bei Gate-Pass: Benchmark-Lauf.
3. Behalten nur bei Verbesserung über der Rausch-Schwelle (~3 %, aus
   Wiederholungsmessungen bestimmt), sonst Revert.
4. Journal-Eintrag immer — auch Fehlschläge. Das Journal ist das Gedächtnis des Loops
   über Sessions hinweg.

### Mutations-Leiter (Aufwand/Risiko aufsteigend)

1. **Knobs**: `nlim`, `nskip`, `looplim`, `throwp`, `ir`, `ctr`, `limitp`,
   dps-Staffelung der Init-Iterationen.
2. **Fork-Code**: Präzisions-Drosselung in frühen Iterationen, redundante
   Neuberechnungen kappen, Loop-Umbau.
3. **Algorithmisch**: Warmstart von niedrigerem dps, Serien-Beschleunigung, bessere
   Startwerte.

### Anti-Gaming-Regeln

- `gate.py` und `reference/` sind vom Loop unantastbar.
- Timing immer auf identischem Workload; Cache-Zustand kontrolliert (kalt/warm
  getrennt).
- Referenzwerte stammen ausschließlich aus der eingefrorenen Original-fatou.gp.

## Nicht-Ziele

- Kein Kernel-Neubau außerhalb von PARI/GP (C/Arb/FLINT) in diesem Vorhaben.
- Keine Änderung der mathematischen Semantik (Kneser-Lösung bleibt die Zielfunktion).
- Keine Anbindung der externen autoresearch-Harness als Voraussetzung.

## Risiken

- **State-Cache**: unklarer Global-Umfang von sexpinit; Mitigation: Anker-Validierung
  nach readbin + Fallback auf Voll-Init.
- **Pipe-Protokoll**: GP-Fehlerausgaben können den Marker-Stream verschmutzen;
  Mitigation: Timeout + Worker-Restart + Rohausgabe im Fehler.
- **Benchmark-Rauschen** (Windows, Hintergrundlast): Wiederholungsmessungen, Median,
  Rausch-Schwelle vor Loop-Start kalibrieren.
- **32→64-bit-Wechsel** kann numerische Details verschieben: Referenzwerte werden nach
  dem Wechsel einmalig neu erzeugt und dann eingefroren.

## Meilensteine

1. **M1**: 64-bit PARI + persistente Session + Benchmark; Messung alt vs. neu.
2. **M2**: State-Cache + Worker-Pool; Messung.
3. **M3**: research/-Struktur (program.md, gate.py, reference/, journal.md,
   fatou_fork.gp); Rausch-Schwelle kalibriert; erste Knob-Experimente.
4. **M4+**: Loop-Betrieb (Knobs → Fork-Code → algorithmisch), auf Wunsch autonom
   via /loop.
