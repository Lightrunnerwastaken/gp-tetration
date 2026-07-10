# Autoresearch-Programm: fatou_fork.gp beschleunigen

## Ziel

Verifizierte korrekte Stellen pro Sekunde maximieren (bench/benchmark.py auf dem
eingefrorenen Workload), ohne das Gate zu verletzen.

## Unantastbar (Anti-Gaming) — Update 2026-07-10

- Nutzer-Direktive (Chat, 2026-07-10): Migrationen/Entscheidungen an
  Gate/Referenzen darf der Loop SELBST treffen, wenn sie nachweislich
  hilfreich sind und die Leistung (Genauigkeit oder Speed) verbessern.
- Die VERIFIKATIONS-Disziplin bleibt unveraendert Pflicht: Jede Referenz-
  Aenderung braucht Dekorrelations-Verifikation (unabhaengige Laeufe mit
  verschiedener Konfiguration muessen ueber die noetige Stellenzahl hinaus
  uebereinstimmen), jede Gate-Schwellen-Aenderung muss wahrheitsbasiert
  kalibriert und im Journal + einem DECISION-Dokument begruendet sein.
  Schwellen duerfen nie unter das sinken, was der aktuelle Fork nachweislich
  echt liefert (kein Selbst-Weichspuelen).
- `bench/workload.py` und `bench/metrics.py` (die Messgroesse selbst)
  bleiben eingefroren; Referenz ist immer die Original-`fatou.gp`.

## Experiment-Protokoll (ein Experiment = eine Mutation)

1. Hypothese im Journal notieren (was wird geaendert, warum sollte es schneller sein).
2. Genau EINE Mutation an `src/fatou_backend/vendor/fatou_fork.gp` ODER an den
   Wrapper-Knobs (nlim/nskip/looplim/dps-Staffelung) vornehmen.
3. `python research/gate.py` -> bei FAIL: revert, Journal-Eintrag, naechstes Experiment.
4. Bei PASS: `python bench/benchmark.py --fatou src/fatou_backend/vendor/fatou_fork.gp
   --label exp-<NNN> --repeat 3` (warm-Median zaehlt; fuer Init-Experimente cold-Median
   mit frischem FATOU_CACHE_DIR).
5. Keep-Regel (erweitert 2026-07-10 auf Nutzer-Direktive: Genauigkeit UND Speed):
   behalten wenn ENTWEDER (a) digits/s um mehr als die Rausch-Schwelle aus
   `research/noise.json` steigt und die Gate-worst-digits nicht sinken, ODER
   (b) die Gate-worst-digits einer Gruppe messbar steigen (>=1 digit, z.B. der
   e|80-Ausreisser) und digits/s nicht um mehr als die Rausch-Schwelle faellt.
   Sonst revert.
6. Journal-Eintrag IMMER (auch Fehlschlaege): Ergebnisdatei, Entscheidung, Learnings.
7. Behaltene Mutation committen (`research: exp-<NNN> <kurzbeschreibung>`).

## Mutations-Leiter (Aufwand/Risiko aufsteigend)

1. Knobs: nlim, nskip, looplim, throwp, ir, ctr, limitp, theta0lim,
   dps-Staffelung der Init-Iterationen.
2. Fork-Code: Praezisions-Drosselung in fruehen sexpinit-Iterationen, redundante
   Neuberechnungen kappen, Schleifen-Restrukturierung.
3. Algorithmisch: Warmstart von niedrigerem dps, Serien-Beschleunigung,
   bessere Startwerte.

## Messdisziplin

- Timing nur ueber bench/benchmark.py (fester Workload, warm/cold getrennt).
- Kaltstart-Vergleiche nur bei kontrolliertem Cache (FATOU_CACHE_DIR auf ein
  frisches Verzeichnis setzen oder state_cache=False).
- Rausch-Schwelle: research/noise.json (aus research/calibrate_noise.py).
- Waehrend Benchmark-Laeufen keine anderen rechenintensiven Prozesse.

## Bekannte Eigenheiten der Engine (Stand Kalibrierung 2026-07-10)

- `e|80` hat einen Ausreisser-Case mit nur ~62.7 Stellen Uebereinstimmung zur
  dps+20-Referenz (auch im Original) -> Gate-Schwelle dort 57.6.
- GP-Identifier duerfen nicht mit Unterstrich beginnen.
- sexpinit-State ist vollstaendig ueber \uv-Variablen + writebin/read
  restaurierbar (Basis des State-Caches).
