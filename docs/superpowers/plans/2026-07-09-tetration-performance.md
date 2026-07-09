# Tetration-Performance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Mehr verifizierte korrekte Stellen pro Sekunde aus dem fatou.gp-Kern: persistente GP-Sessions, 64-bit PARI, sexpinit-State-Cache, Worker-Pool, plus Benchmark/Gate/research-Struktur für den Autoresearch-Loop (Meilensteine M1–M3 der Spec `docs/superpowers/specs/2026-07-09-tetration-performance-design.md`).

**Architecture:** Der bestehende One-Shot-Pfad (`gp.exe` pro Batch neu starten) bleibt als Fallback erhalten; ein neuer `FatouGPWorker` hält einen langlebigen `gp.exe`-Prozess pro (Basis, dps) mit markerbasiertem Pipe-Protokoll. Ein State-Cache (`writebin`/`read`) macht `sexpinit` pro (Basis, dps, Knobs, fatou-Hash) einmalig. Benchmark und Gate messen alles gegen eingefrorene Referenzwerte der Original-fatou.gp.

**Tech Stack:** Python ≥3.11, mpmath (einzige Dependency), PARI/GP 2.17.3 (`gp.exe`), unittest (Repo-Konvention, kein pytest).

## Global Constraints

- Repo: `C:\Users\janis\Documents\Tetration\gp-tetration`, Branch `tetration-performance`.
- Python ≥3.11; einzige Runtime-Dependency bleibt `mpmath>=1.3` (keine neuen Dependencies).
- Testrunner ist **unittest** mit `sys.path.insert(0, …/src)`-Muster wie in `tests/test_gp_backend.py`; Schnell-Suite: `python -m unittest discover -s tests -p "test_*.py" -v` (ohne `FATOU_BACKEND_RUN_SLOW`) muss nach jedem Task grün sein.
- gp-spawnende Tests benutzen kleine Parameter (`dps=50, nlim=20, nskip=4, looplim=30`), damit die Suite unter ~2 Minuten bleibt.
- GP-Skript-Pfade immer mit `Path(...).as_posix()` (Forward-Slashes) in GP-Strings einbetten.
- Die öffentliche API (`FatouGP`, `FatouGPSession`, CLI-Flags) darf nicht brechen; `mixed_phase` und `cli` laufen unverändert weiter.
- Ab Task 9 sind `research/gate.py` und `research/reference/` eingefroren (Anti-Gaming): kein späterer Task und kein Loop-Experiment darf sie ändern.
- Commits: Conventional-Prefix (feat/test/docs/chore), Message-Text Deutsch, jede Commit-Message endet mit `Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>`.
- Code und Code-Kommentare Englisch (Repo-Konvention); `research/program.md` und `research/journal.md` Deutsch.

**Vorab verifizierte Fakten** (am 2026-07-09 gegen das installierte `gp.exe` 2.17.3 32-bit geprüft — nicht erneut in Frage stellen):

1. Die bestehende Zeile `\\p {dps}` setzt die realprecision korrekt (Ausgabe „realprecision = 86 significant digits (80 digits displayed)" bei dps 80). Kein Bug; das Refactoring auf `default(realprecision, N);` ist nur Klarheit, keine Fehlerkorrektur.
2. `writebin("f", [a,b,…])` schreibt PARI-Objekte binär; zurücklesen mit `vv = read("f")` (es gibt **kein** `readbin`). Präzision bleibt erhalten (86 Digits nach Roundtrip verifiziert).
3. Der Metabefehl `\uv` listet alle User-Variablen im Format `name =\n  <wert>` — Variablennamen matchen auf Regex `^([A-Za-z_][A-Za-z0-9_]*) =$`, Wertzeilen sind eingerückt.
4. 32-bit-PARI default realprecision rundet auf Limb-Vielfache auf (80 → 86 Digits).

---

### Task 1: `_init_lines`-Refactor (gemeinsame Init-Sequenz für One-Shot und Worker)

**Files:**
- Modify: `src/fatou_backend/gp_backend.py` (Methode `_run_initialized`, neue Methode `_init_lines`)
- Test: `tests/test_init_lines.py` (neu)

**Interfaces:**
- Produces: `FatouGP._init_lines(base: GPValue) -> list[str]` — die GP-Zeilen, die einen frischen gp-Prozess bis inkl. `sexpinit` initialisieren. Task 5 (Worker) und Task 7 (Cache) bauen darauf auf.

- [ ] **Step 1: Failing Test schreiben**

`tests/test_init_lines.py`:

```python
from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend import FatouGP, find_default_fatou_gp, find_default_gp_exe


class InitLinesTests(unittest.TestCase):
    def setUp(self) -> None:
        self.gp = FatouGP(
            gp_exe=find_default_gp_exe(),
            fatou_gp=find_default_fatou_gp(),
            dps=50,
            nlim=20,
            nskip=4,
            looplim=30,
        )

    def test_init_lines_order_and_content(self) -> None:
        lines = self.gp._init_lines("e")
        self.assertEqual(lines[0], "default(realprecision, 50);")
        self.assertTrue(lines[1].startswith('read("'))
        self.assertTrue(lines[1].endswith('");'))
        self.assertIn("fatou.gp", lines[1])
        self.assertNotIn("\\", lines[1])  # as_posix, keine Backslashes
        self.assertEqual(lines[2], "quietmode=1;")
        self.assertEqual(lines[3], "sexpinit(exp(1),20,4,30);")
        self.assertEqual(len(lines), 4)

    def test_init_lines_complex_base(self) -> None:
        lines = self.gp._init_lines("1+I")
        self.assertEqual(lines[3], "sexpinit(1+I,20,4,30);")


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Test laufen lassen, Fehlschlag verifizieren**

Run: `cd "C:\Users\janis\Documents\Tetration\gp-tetration" && python -m unittest tests.test_init_lines -v`
Expected: FAIL/ERROR mit `AttributeError: 'FatouGP' object has no attribute '_init_lines'`

- [ ] **Step 3: Implementieren**

In `src/fatou_backend/gp_backend.py`, direkt über `_run_initialized` einfügen:

```python
    def _init_lines(self, base: GPValue) -> list[str]:
        fatou_path = Path(self.fatou_gp).as_posix()
        return [
            f"default(realprecision, {self.dps});",
            f'read("{fatou_path}");',
            f"quietmode={self.quietmode};",
            f"sexpinit({self._base_expr(base)},{self.nlim},{self.nskip},{self.looplim});",
        ]
```

Und in `_run_initialized` den Zeilen-Aufbau ersetzen — aus:

```python
        fatou_path = self.fatou_gp.as_posix()
        lines = [
            f"\\\\p {self.dps}",
            f'read("{fatou_path}")',
            f"quietmode={self.quietmode};",
            f"sexpinit({self._base_expr(base)},{self.nlim},{self.nskip},{self.looplim});",
            'print("__BEGIN_RESULTS__")',
        ]
```

wird:

```python
        lines = self._init_lines(base) + ['print("__BEGIN_RESULTS__")']
```

- [ ] **Step 4: Alle Tests laufen lassen**

Run: `cd "C:\Users\janis\Documents\Tetration\gp-tetration" && python -m unittest discover -s tests -p "test_*.py" -v`
Expected: alle Tests PASS (12 bestehende + 2 neue). Die gp-Tests beweisen, dass `default(realprecision, …)` + Semikolon hinter `read(…)` das Verhalten nicht ändern.

- [ ] **Step 5: Commit**

```bash
git add src/fatou_backend/gp_backend.py tests/test_init_lines.py
git commit -m "refactor: Init-Sequenz in _init_lines extrahiert (Basis fuer Worker/Cache)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 2: 64-bit PARI bevorzugen und installieren

**Files:**
- Modify: `src/fatou_backend/gp_backend.py` (`find_default_gp_exe`, neues Modul-Konstant `GP_EXE_CANDIDATES`)
- Test: `tests/test_init_lines.py` (erweitern)

**Interfaces:**
- Produces: `GP_EXE_CANDIDATES: tuple[Path, ...]` (Modul-Level, Pari64 zuerst). `find_default_gp_exe()` Signatur unverändert.

- [ ] **Step 1: Failing Test schreiben**

In `tests/test_init_lines.py` ergänzen (Import erweitern um `from fatou_backend.gp_backend import GP_EXE_CANDIDATES`):

```python
class GpExeCandidateTests(unittest.TestCase):
    def test_pari64_is_preferred(self) -> None:
        self.assertIn("Pari64", str(GP_EXE_CANDIDATES[0]))
        self.assertIn("Pari32", str(GP_EXE_CANDIDATES[1]))
```

- [ ] **Step 2: Test laufen lassen, Fehlschlag verifizieren**

Run: `python -m unittest tests.test_init_lines -v`
Expected: ERROR `ImportError: cannot import name 'GP_EXE_CANDIDATES'`

- [ ] **Step 3: Implementieren**

In `src/fatou_backend/gp_backend.py` — aus dem Funktionskörper von `find_default_gp_exe` die Kandidaten auf Modul-Ebene ziehen und Reihenfolge drehen:

```python
GP_EXE_CANDIDATES = (
    Path(r"C:\Program Files\Pari64-2-17-3\gp.exe"),
    Path(r"C:\Program Files (x86)\Pari32-2-17-3\gp.exe"),
)


def find_default_gp_exe() -> Path:
    env = os.getenv("FATOU_GP_EXE")
    if env:
        path = Path(env)
        if path.exists():
            return path
    path = _first_existing(GP_EXE_CANDIDATES)
    if path is None:
        raise FileNotFoundError("Could not locate gp.exe. Set FATOU_GP_EXE.")
    return path
```

- [ ] **Step 4: Tests laufen lassen**

Run: `python -m unittest tests.test_init_lines -v`
Expected: PASS

- [ ] **Step 5: Pari64 installieren (Umgebungsschritt, PowerShell)**

```powershell
Invoke-WebRequest "https://pari.math.u-bordeaux.fr/pub/pari/windows/Pari64-2-17-3.exe" -OutFile "$env:TEMP\Pari64-setup.exe"
Start-Process "$env:TEMP\Pari64-setup.exe" -ArgumentList "/S" -Wait
& "C:\Program Files\Pari64-2-17-3\gp.exe" --version
```

Expected: Versionsbanner mit `Version 2.17.3` und `amd64`.
**Fallback:** Wenn der Download 404 liefert, im Verzeichnislisting `https://pari.math.u-bordeaux.fr/pub/pari/windows/` den aktuellsten `Pari64-2-17-*.exe` nehmen; wenn dessen Installationspfad vom Kandidaten abweicht, den tatsächlichen Pfad zusätzlich VOR den bestehenden Einträgen in `GP_EXE_CANDIDATES` eintragen (Test in Step 1 ggf. auf den neuen Pfadnamen anpassen). Wenn gar kein 64-bit-Download möglich ist: Schritt überspringen, Notiz in der Commit-Message — der Code-Teil dieses Tasks ist trotzdem korrekt (Fallback auf Pari32 funktioniert).

- [ ] **Step 6: Volle Suite gegen 64-bit laufen lassen**

Run: `python -m unittest discover -s tests -p "test_*.py" -v`
Expected: alle PASS. `find_default_gp_exe()` findet jetzt Pari64 (prüfen: `python -c "import sys; sys.path.insert(0,'src'); from fatou_backend import find_default_gp_exe; print(find_default_gp_exe())"` → Pfad enthält `Pari64`).

- [ ] **Step 7: Commit**

```bash
git add src/fatou_backend/gp_backend.py tests/test_init_lines.py
git commit -m "feat: 64-bit PARI bevorzugt (Pari64 vor Pari32 in Kandidatenliste)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 3: Workload-Definition und Referenzwerte

**Files:**
- Create: `bench/__init__.py` (leer), `bench/workload.py`, `bench/make_reference.py`
- Create (generiert): `research/reference/values.json`
- Test: `tests/test_workload.py`

**Interfaces:**
- Produces:
  - `bench.workload.Case` — frozen dataclass: `kind: str` ("sexp"|"slog"), `base: str`, `arg: str`, `dps: int`, `tier: str` ("throughput"|"precision"|"anchor").
  - `bench.workload.all_cases(deep: bool = False) -> list[Case]`
  - `bench.workload.case_id(c: Case) -> str` — eindeutiger Schlüssel, Format `f"{c.kind}|{c.base}|{c.arg}|{c.dps}"`.
  - `research/reference/values.json` — Schema:
    ```json
    {"meta": {"generated": "...", "gp_exe": "...", "fatou_sha256": "...", "ref_dps_margin": 20},
     "values": {"sexp|e|0.5|80": {"real": "...", "imag": "..."}, "...": {}}}
    ```
- Consumes: `FatouGP` aus Task 1/2.

- [ ] **Step 1: Failing Test schreiben**

`tests/test_workload.py`:

```python
from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from bench.workload import Case, all_cases, case_id


class WorkloadTests(unittest.TestCase):
    def test_default_cases_shape(self) -> None:
        cases = all_cases()
        tiers = {c.tier for c in cases}
        self.assertEqual(tiers, {"throughput", "precision", "anchor"})
        throughput = [c for c in cases if c.tier == "throughput"]
        self.assertEqual({c.base for c in throughput}, {"e", "2", "10"})
        self.assertTrue(all(c.dps == 80 for c in throughput))
        # 16 sexp- + 16 slog-Argumente pro Basis
        self.assertEqual(len(throughput), 3 * 32)
        precision = [c for c in cases if c.tier == "precision"]
        self.assertEqual({(c.base, c.dps) for c in precision}, {("e", 200), ("2", 200)})
        anchors = [c for c in cases if c.tier == "anchor"]
        self.assertEqual({c.base for c in anchors}, {"1+I", "0.8+0.4*I", "2+I"})

    def test_deep_adds_500_and_1000(self) -> None:
        deep = [c for c in all_cases(deep=True) if c.tier == "precision"]
        self.assertEqual({c.dps for c in deep}, {200, 500, 1000})

    def test_case_id_roundtrip_unique(self) -> None:
        cases = all_cases(deep=True)
        ids = [case_id(c) for c in cases]
        self.assertEqual(len(ids), len(set(ids)))


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Test laufen lassen, Fehlschlag verifizieren**

Run: `python -m unittest tests.test_workload -v`
Expected: ERROR `ModuleNotFoundError: No module named 'bench'`

- [ ] **Step 3: Implementieren**

`bench/__init__.py`: leere Datei.

`bench/workload.py`:

```python
"""Fixed benchmark/gate workload. Frozen once reference values exist."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Case:
    kind: str   # "sexp" | "slog"
    base: str   # GP base expression, e.g. "e", "2", "1+I"
    arg: str    # GP-parsable argument literal
    dps: int
    tier: str   # "throughput" | "precision" | "anchor"


THROUGHPUT_BASES = ("e", "2", "10")
THROUGHPUT_DPS = 80
# 16 sexp arguments in [-0.4375, 1.4375]
SEXP_ARGS = tuple(str(-7 / 16 + j / 8) for j in range(16))
# 16 slog arguments in [1.25, 5.0]
SLOG_ARGS = tuple(str(1.25 + j / 4) for j in range(16))
PRECISION_CASES = (("e", "0.5"), ("2", "0.5"))
PRECISION_DPS = (200,)
PRECISION_DPS_DEEP = (200, 500, 1000)
ANCHOR_BASES = ("1+I", "0.8+0.4*I", "2+I")
ANCHOR_DPS = 80


def all_cases(deep: bool = False) -> list[Case]:
    cases: list[Case] = []
    for base in THROUGHPUT_BASES:
        for arg in SEXP_ARGS:
            cases.append(Case("sexp", base, arg, THROUGHPUT_DPS, "throughput"))
        for arg in SLOG_ARGS:
            cases.append(Case("slog", base, arg, THROUGHPUT_DPS, "throughput"))
    for base, arg in PRECISION_CASES:
        for dps in (PRECISION_DPS_DEEP if deep else PRECISION_DPS):
            cases.append(Case("sexp", base, arg, dps, "precision"))
    for base in ANCHOR_BASES:
        cases.append(Case("sexp", base, "0.5", ANCHOR_DPS, "anchor"))
    return cases


def case_id(c: Case) -> str:
    return f"{c.kind}|{c.base}|{c.arg}|{c.dps}"
```

`bench/make_reference.py`:

```python
"""Generate frozen reference values with the ORIGINAL fatou.gp at dps+margin.

Run once after the 64-bit switch; afterwards research/reference/ is frozen.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend import FatouGP, find_default_fatou_gp, find_default_gp_exe

from workload import all_cases, case_id

REF_MARGIN = 20
REFERENCE_PATH = Path(__file__).resolve().parents[1] / "research" / "reference" / "values.json"


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate frozen reference values")
    parser.add_argument("--deep", action="store_true", help="include dps 500/1000 tiers")
    parser.add_argument("--out", default=str(REFERENCE_PATH))
    args = parser.parse_args()

    fatou = find_default_fatou_gp()
    gp_exe = find_default_gp_exe()
    values: dict[str, dict[str, str]] = {}
    # group by (base, dps) so each sexpinit happens once per group
    groups: dict[tuple[str, int], list] = defaultdict(list)
    for c in all_cases(deep=args.deep):
        groups[(c.base, c.dps)].append(c)

    for (base, dps), cases in sorted(groups.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        ref_dps = dps + REF_MARGIN
        mp.mp.dps = ref_dps + 30
        gp = FatouGP(gp_exe=gp_exe, fatou_gp=fatou, dps=ref_dps,
                     looplim=max(35, ref_dps - 20))
        expressions = [f"{c.kind}({c.arg})" for c in cases]
        print(f"[make_reference] base={base} dps={dps} (ref_dps={ref_dps}) "
              f"{len(expressions)} exprs ...", flush=True)
        results = gp.eval_batch(base, expressions)
        for c, val in zip(cases, results):
            values[case_id(c)] = {
                "real": mp.nstr(mp.re(val), ref_dps, min_fixed=0, max_fixed=0),
                "imag": mp.nstr(mp.im(val), ref_dps, min_fixed=0, max_fixed=0),
            }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "meta": {
            "generated": datetime.date.today().isoformat(),
            "gp_exe": str(gp_exe),
            "fatou_sha256": hashlib.sha256(Path(fatou).read_bytes()).hexdigest(),
            "ref_dps_margin": REF_MARGIN,
        },
        "values": values,
    }
    out.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"[make_reference] wrote {len(values)} values to {out}")


if __name__ == "__main__":
    main()
```

Hinweis: `from workload import …` funktioniert, weil das Skript direkt (`python bench/make_reference.py`) läuft und sein eigenes Verzeichnis auf `sys.path` liegt.

- [ ] **Step 4: Tests laufen lassen**

Run: `python -m unittest tests.test_workload -v`
Expected: PASS (3 Tests)

- [ ] **Step 5: Referenzwerte generieren (ohne --deep)**

Run: `cd "C:\Users\janis\Documents\Tetration\gp-tetration" && python bench/make_reference.py`
Expected: Fortschrittszeilen pro (base, dps)-Gruppe; am Ende `wrote 101 values to ...` (96 throughput [3 Basen × 32] + 2 precision + 3 anchor = 101 = `len(all_cases())`).
Laufzeit: Minuten-Bereich; die dps-220-Gruppen sind die langsamsten. `--deep` (500/1000) wird ERST nach M2 nachgezogen (separater Lauf, gleiche Datei wird ergänzt — dazu Task-10-Hinweis).

- [ ] **Step 6: Stichprobe gegen bekannte Werte prüfen**

Run: `python -c "import json; v=json.load(open('research/reference/values.json'))['values']; print(v['sexp|e|0.5|80']['real'][:22])"`
Expected: beginnt mit `1.64635423375119458097` (bekannter Wert aus `tests/test_gp_backend.py`).

- [ ] **Step 7: Commit**

```bash
git add bench/__init__.py bench/workload.py bench/make_reference.py research/reference/values.json tests/test_workload.py
git commit -m "feat: Benchmark-Workload + eingefrorene Referenzwerte (Original-fatou.gp, dps+20)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 4: Digit-Metrik, Benchmark-CLI und Baseline-Messung

**Files:**
- Create: `bench/metrics.py`, `bench/benchmark.py`, `bench/results/` (Verzeichnis)
- Test: `tests/test_metrics.py`

**Interfaces:**
- Produces:
  - `bench.metrics.correct_digits(value: mp.mpc, reference: mp.mpc, cap: int) -> float` — `min(cap, -log10(|v-ref| / max(1,|ref|)))`, `cap` bei exakter Übereinstimmung.
  - `bench.metrics.load_reference(path: Path) -> dict[str, mp.mpc]` — case_id → mpc (setzt `mp.mp.dps` NICHT selbst; Aufrufer muss dps hoch genug halten).
  - CLI: `python bench/benchmark.py --mode all|throughput|precision --repeat N --deep --fatou PATH --label TEXT --json PATH` — schreibt Ergebnis-JSON, druckt Zusammenfassung.
  - Ergebnis-JSON-Schema (von Task 6/8/10 konsumiert):
    ```json
    {"meta": {"label": "...", "date": "...", "gp_exe": "...", "fatou": "...", "repeat": 1},
     "runs": [{"group": "throughput|e|80", "mode": "cold", "seconds": 1.23,
               "values": 32, "min_correct_digits": 71.2, "digits_per_second": 1850.0}]}
    ```
- Consumes: Task 3 (`workload`, `values.json`), `FatouGP`.

- [ ] **Step 1: Failing Test für die Metrik schreiben**

`tests/test_metrics.py`:

```python
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from bench.metrics import correct_digits


class MetricsTests(unittest.TestCase):
    def setUp(self) -> None:
        mp.mp.dps = 120

    def test_exact_match_returns_cap(self) -> None:
        v = mp.mpc("1.5", "0.25")
        self.assertEqual(correct_digits(v, v, cap=80), 80.0)

    def test_known_relative_error(self) -> None:
        ref = mp.mpc("2.0", "0")
        val = ref + mp.mpf("2e-10")  # rel err 1e-10
        d = correct_digits(val, ref, cap=80)
        self.assertAlmostEqual(d, 10.0, places=3)

    def test_small_reference_uses_absolute_scale(self) -> None:
        ref = mp.mpc("1e-30", "0")   # |ref| < 1 -> scale = 1
        val = ref + mp.mpf("1e-12")
        d = correct_digits(val, ref, cap=80)
        self.assertAlmostEqual(d, 12.0, places=3)

    def test_cap_is_upper_bound(self) -> None:
        ref = mp.mpc("2.0", "0")
        val = ref + mp.mpf("1e-300")
        self.assertEqual(correct_digits(val, ref, cap=50), 50.0)


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Test laufen lassen, Fehlschlag verifizieren**

Run: `python -m unittest tests.test_metrics -v`
Expected: ERROR `ModuleNotFoundError: No module named 'bench.metrics'`

- [ ] **Step 3: Metrik implementieren**

`bench/metrics.py`:

```python
"""Digit-accuracy metric and reference loading. Pure functions, no GP calls."""
from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp


def correct_digits(value: mp.mpc, reference: mp.mpc, cap: int) -> float:
    err = abs(mp.mpc(value) - mp.mpc(reference))
    if err == 0:
        return float(cap)
    scale = max(mp.mpf(1), abs(mp.mpc(reference)))
    digits = -mp.log10(err / scale)
    return float(min(mp.mpf(cap), digits))


def load_reference(path: Path) -> dict[str, mp.mpc]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return {
        key: mp.mpc(mp.mpf(entry["real"]), mp.mpf(entry["imag"]))
        for key, entry in payload["values"].items()
    }
```

- [ ] **Step 4: Metrik-Tests laufen lassen**

Run: `python -m unittest tests.test_metrics -v`
Expected: PASS (4 Tests)

- [ ] **Step 5: Benchmark-CLI implementieren**

`bench/benchmark.py`:

```python
"""Benchmark: verified correct digits per second on the frozen workload.

cold  = fresh FatouGP per group, timing includes process start + sexpinit
warm  = same FatouGP object, second evaluation of the group (Task 6 makes
        this meaningful via persistent workers; on the one-shot path warm
        equals cold and is reported anyway for comparability)
"""
from __future__ import annotations

import argparse
import datetime
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend import FatouGP, find_default_fatou_gp, find_default_gp_exe

from metrics import correct_digits, load_reference
from workload import all_cases, case_id

REFERENCE_PATH = Path(__file__).resolve().parents[1] / "research" / "reference" / "values.json"
RESULTS_DIR = Path(__file__).resolve().parent / "results"


def run_group(gp: FatouGP, base: str, cases: list, reference: dict) -> dict:
    expressions = [f"{c.kind}({c.arg})" for c in cases]
    start = time.perf_counter()
    results = gp.eval_batch(base, expressions)
    elapsed = time.perf_counter() - start
    dps = cases[0].dps
    digits = [correct_digits(v, reference[case_id(c)], cap=dps - 10)
              for c, v in zip(cases, results)]
    min_digits = min(digits)
    return {
        "seconds": elapsed,
        "values": len(cases),
        "min_correct_digits": min_digits,
        "digits_per_second": sum(digits) / elapsed,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="fatou_backend benchmark")
    parser.add_argument("--mode", choices=["all", "throughput", "precision"], default="all")
    parser.add_argument("--repeat", type=int, default=1)
    parser.add_argument("--deep", action="store_true")
    parser.add_argument("--fatou", default=None, help="fatou.gp variant to benchmark")
    parser.add_argument("--label", default="run")
    parser.add_argument("--json", default=None)
    args = parser.parse_args()

    fatou = Path(args.fatou) if args.fatou else find_default_fatou_gp()
    gp_exe = find_default_gp_exe()
    reference_raw = load_reference(REFERENCE_PATH)

    tiers = {"all": {"throughput", "precision"}, "throughput": {"throughput"},
             "precision": {"precision"}}[args.mode]
    groups: dict[tuple[str, int], list] = defaultdict(list)
    for c in all_cases(deep=args.deep):
        if c.tier in tiers:
            groups[(c.base, c.dps)].append(c)

    runs = []
    for (base, dps), cases in sorted(groups.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        mp.mp.dps = dps + 40
        reference = load_reference(REFERENCE_PATH)  # reparse at sufficient dps
        for rep in range(args.repeat):
            gp = FatouGP(gp_exe=gp_exe, fatou_gp=fatou, dps=dps,
                         looplim=max(35, dps - 20))
            tier = cases[0].tier
            cold = run_group(gp, base, cases, reference)
            cold.update({"group": f"{tier}|{base}|{dps}", "mode": "cold", "rep": rep})
            runs.append(cold)
            warm = run_group(gp, base, cases, reference)
            warm.update({"group": f"{tier}|{base}|{dps}", "mode": "warm", "rep": rep})
            runs.append(warm)
            close = getattr(gp, "close", None)
            if close is not None:
                close()
            print(f"{tier}|{base}|{dps} rep {rep}: "
                  f"cold {cold['seconds']:.2f}s ({cold['digits_per_second']:.0f} d/s), "
                  f"warm {warm['seconds']:.2f}s ({warm['digits_per_second']:.0f} d/s), "
                  f"min digits {min(cold['min_correct_digits'], warm['min_correct_digits']):.1f}",
                  flush=True)

    payload = {
        "meta": {"label": args.label, "date": datetime.datetime.now().isoformat(timespec="seconds"),
                 "gp_exe": str(gp_exe), "fatou": str(fatou), "repeat": args.repeat},
        "runs": runs,
    }
    out = Path(args.json) if args.json else RESULTS_DIR / f"{datetime.date.today().isoformat()}-{args.label}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"[benchmark] wrote {out}")


if __name__ == "__main__":
    main()
```

Anmerkung: `load_reference` wird pro Gruppe nach dem Setzen von `mp.mp.dps` erneut aufgerufen, damit die Referenz-mpc-Werte mit genug Stellen geparst sind; `reference_raw` dient nur dem Fail-Fast beim Start (Datei fehlt → sofortiger Fehler). `gp.close()` existiert erst ab Task 6, daher der `getattr`-Guard.

- [ ] **Step 6: Baseline messen (One-Shot-Welt, vor Worker)**

Run: `python bench/benchmark.py --mode all --label baseline-oneshot`
Expected: pro Gruppe eine Zeile; cold ≈ warm (One-Shot-Pfad hat keinen Warm-Vorteil); `min digits` ≥ 60 bei dps 80 (d. h. Referenz-Übereinstimmung weit über Rauschen). JSON landet in `bench/results/2026-07-09-baseline-oneshot.json` (Datum entsprechend).

- [ ] **Step 7: Schnelle Suite laufen lassen**

Run: `python -m unittest discover -s tests -p "test_*.py" -v`
Expected: alle PASS

- [ ] **Step 8: Commit**

```bash
git add bench/metrics.py bench/benchmark.py bench/results/ tests/test_metrics.py
git commit -m "feat: Digit-Metrik + Benchmark-CLI + One-Shot-Baseline-Messung

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 5: FatouGPWorker — persistenter GP-Prozess

**Files:**
- Create: `src/fatou_backend/worker.py`
- Test: `tests/test_worker.py`

**Interfaces:**
- Consumes: `FatouGP._init_lines(base)` (Task 1) — der Worker selbst ist aber von `FatouGP` entkoppelt und nimmt rohe Zeilen.
- Produces:
  - `class WorkerDied(RuntimeError)` — Prozess tot/Timeout; Aufrufer darf genau einmal respawnen.
  - `class FatouGPWorker`:
    - `__init__(self, gp_exe: Path | str, init_lines: list[str], init_timeout: float = 3600.0, eval_timeout: float = 600.0)` — spawnt und initialisiert blockierend.
    - `eval(self, expressions: list[str]) -> list[mp.mpc]`
    - `eval_raw(self, lines: list[str], sentinel: str, timeout: float) -> list[str]` — sendet rohe GP-Zeilen, sammelt Output bis `sentinel`-Zeile (für `\uv`-Discovery in Task 7).
    - `close(self) -> None`, `alive: bool` (Property), Context-Manager (`__enter__`/`__exit__`), `__del__`-Fallback ruft `close()`.

- [ ] **Step 1: Failing Tests schreiben**

`tests/test_worker.py`:

```python
from __future__ import annotations

import sys
import unittest
from pathlib import Path

from mpmath import mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend import FatouGP, find_default_fatou_gp, find_default_gp_exe
from fatou_backend.worker import FatouGPWorker, WorkerDied


def make_gp() -> FatouGP:
    return FatouGP(gp_exe=find_default_gp_exe(), fatou_gp=find_default_fatou_gp(),
                   dps=50, nlim=20, nskip=4, looplim=30)


class WorkerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        mp.dps = 80
        cls.gp = make_gp()
        cls.worker = FatouGPWorker(cls.gp.gp_exe, cls.gp._init_lines("e"))

    @classmethod
    def tearDownClass(cls) -> None:
        cls.worker.close()

    def test_eval_matches_oneshot(self) -> None:
        got = self.worker.eval(["sexp(0.5)"])[0]
        want = self.gp.sexp("e", mp.mpf("0.5"))
        self.assertLess(abs(got - want), mp.mpf("1e-30"))

    def test_second_eval_reuses_process(self) -> None:
        pid_before = self.worker._proc.pid
        values = self.worker.eval(["sexp(0.25)", "slog(2.0)"])
        self.assertEqual(len(values), 2)
        self.assertEqual(self.worker._proc.pid, pid_before)
        self.assertTrue(self.worker.alive)

    def test_gp_error_raises_with_output(self) -> None:
        with self.assertRaises(RuntimeError) as ctx:
            self.worker.eval(["undefined_function_xyz(1)"])
        self.assertIn("***", str(ctx.exception))
        # worker was killed on error; next use must fail fast
        self.assertFalse(self.worker.alive)
        with self.assertRaises(WorkerDied):
            self.worker.eval(["sexp(0.5)"])
        # respawn a shared worker for remaining tests in this class
        type(self).worker = FatouGPWorker(self.gp.gp_exe, self.gp._init_lines("e"))

    def test_timeout_kills_worker(self) -> None:
        w = FatouGPWorker(self.gp.gp_exe, self.gp._init_lines("e"), eval_timeout=3.0)
        with self.assertRaises(WorkerDied):
            w.eval(["while(1,)"])
        self.assertFalse(w.alive)
        w.close()

    def test_context_manager_closes(self) -> None:
        with FatouGPWorker(self.gp.gp_exe, self.gp._init_lines("e")) as w:
            w.eval(["sexp(0.1)"])
            self.assertTrue(w.alive)
        self.assertFalse(w.alive)


if __name__ == "__main__":
    unittest.main()
```

Hinweis zur Testreihenfolge: unittest läuft alphabetisch; `test_gp_error_raises_with_output` respawnt den Klassen-Worker selbst, damit nachfolgende Tests einen lebenden Worker sehen. (Reihenfolge: context_manager → eval_matches → gp_error → second_eval → timeout; `second_eval` läuft nach `gp_error` und braucht den Respawn.)

- [ ] **Step 2: Tests laufen lassen, Fehlschlag verifizieren**

Run: `python -m unittest tests.test_worker -v`
Expected: ERROR `ModuleNotFoundError: No module named 'fatou_backend.worker'`

- [ ] **Step 3: Implementieren**

`src/fatou_backend/worker.py`:

```python
"""Long-lived gp.exe process: init once (sexpinit), then stream evaluations."""
from __future__ import annotations

import queue
import subprocess
import threading
import time
from pathlib import Path

import mpmath as mp


class WorkerDied(RuntimeError):
    """The gp process is gone (crash, EOF, or timeout kill)."""


class _ReaderThread(threading.Thread):
    def __init__(self, stream, out_queue: "queue.Queue[str | None]") -> None:
        super().__init__(daemon=True)
        self._stream = stream
        self._queue = out_queue

    def run(self) -> None:
        for line in iter(self._stream.readline, ""):
            self._queue.put(line.rstrip("\r\n"))
        self._queue.put(None)  # EOF sentinel


def _parse_gp_scalar(text: str) -> mp.mpf:
    return mp.mpf(text.strip().replace(" ", ""))


class FatouGPWorker:
    def __init__(self, gp_exe: Path | str, init_lines: list[str],
                 init_timeout: float = 3600.0, eval_timeout: float = 600.0) -> None:
        self.eval_timeout = eval_timeout
        self._tag = 0
        self._alive = False
        self._proc = subprocess.Popen(
            [str(gp_exe), "-q"],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1,
        )
        self._queue: "queue.Queue[str | None]" = queue.Queue()
        self._reader = _ReaderThread(self._proc.stdout, self._queue)
        self._reader.start()
        self._alive = True
        self._send(init_lines + ['print("__READY__")'])
        self._read_until("__READY__", init_timeout)

    # -- lifecycle ---------------------------------------------------------

    @property
    def alive(self) -> bool:
        return self._alive and self._proc.poll() is None

    def close(self) -> None:
        if getattr(self, "_proc", None) is None:
            return
        self._alive = False
        try:
            if self._proc.poll() is None:
                self._proc.kill()
                self._proc.wait(timeout=10)
        except Exception:
            pass

    def __enter__(self) -> "FatouGPWorker":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def __del__(self) -> None:
        self.close()

    # -- protocol ----------------------------------------------------------

    def _send(self, lines: list[str]) -> None:
        if not self.alive:
            raise WorkerDied("gp process is not running")
        try:
            self._proc.stdin.write("\n".join(lines) + "\n")
            self._proc.stdin.flush()
        except OSError as exc:
            self.close()
            raise WorkerDied(f"failed to write to gp: {exc}") from exc

    def _read_until(self, sentinel: str, timeout: float) -> list[str]:
        deadline = time.perf_counter() + timeout
        collected: list[str] = []
        errors: list[str] = []
        while True:
            remaining = deadline - time.perf_counter()
            if remaining <= 0:
                self.close()
                raise WorkerDied(
                    f"timeout after {timeout}s waiting for {sentinel!r}. "
                    f"Last output:\n" + "\n".join(collected[-20:]))
            try:
                line = self._queue.get(timeout=min(remaining, 1.0))
            except queue.Empty:
                continue
            if line is None:
                self.close()
                raise WorkerDied(
                    "gp process ended unexpectedly. Last output:\n"
                    + "\n".join(collected[-20:]))
            if line.strip() == sentinel:
                if errors:
                    self.close()
                    raise RuntimeError("PARI/GP error:\n" + "\n".join(errors))
                return collected
            collected.append(line)
            if "***" in line:
                errors.append(line)

    def eval_raw(self, lines: list[str], sentinel: str, timeout: float) -> list[str]:
        self._send(lines + [f'print("{sentinel}")'])
        return self._read_until(sentinel, timeout)

    def eval(self, expressions: list[str]) -> list[mp.mpc]:
        self._tag += 1
        sentinel = f"__END__{self._tag}"
        lines: list[str] = []
        for idx, expr in enumerate(expressions):
            lines.append(f"vv = ({expr});")
            lines.append(f'print("__RES__{idx}")')
            lines.append("print(real(vv));")
            lines.append("print(imag(vv));")
        output = self.eval_raw(lines, sentinel, self.eval_timeout)

        parsed: list[mp.mpc] = []
        cursor = 0
        while cursor < len(output):
            marker = output[cursor].strip()
            if not marker.startswith("__RES__"):
                cursor += 1
                continue
            if cursor + 2 > len(output) - 1:
                break
            real_text = output[cursor + 1]
            imag_text = output[cursor + 2]
            parsed.append(mp.mpc(_parse_gp_scalar(real_text), _parse_gp_scalar(imag_text)))
            cursor += 3
        if len(parsed) != len(expressions):
            self.close()
            raise RuntimeError(
                f"Expected {len(expressions)} results, got {len(parsed)}.\n"
                + "\n".join(output))
        return parsed
```

(Der Parser hat identische Semantik wie der bestehende in `_run_initialized`: Marker suchen, zwei Folgezeilen als real/imag parsen, unvollständige Tripel am Ende brechen ab.)

- [ ] **Step 4: Tests laufen lassen**

Run: `python -m unittest tests.test_worker -v`
Expected: PASS (5 Tests). Der Timeout-Test braucht ~3–4 s, die Klassen-Inits je ~2 s.

- [ ] **Step 5: Volle Suite**

Run: `python -m unittest discover -s tests -p "test_*.py" -v`
Expected: alle PASS

- [ ] **Step 6: Commit**

```bash
git add src/fatou_backend/worker.py tests/test_worker.py
git commit -m "feat: FatouGPWorker - persistenter gp-Prozess mit Pipe-Protokoll, Timeout und Restart-Semantik

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 6: FatouGP-Integration (persistent als Default) + M1-Messung

**Files:**
- Modify: `src/fatou_backend/gp_backend.py` (`FatouGP`: neue Felder, `_worker_for`, Routing in `_run_initialized`, `close`, Context-Manager)
- Modify: `src/fatou_backend/__init__.py` (Export `WorkerDied` optional — nein: unverändert lassen, YAGNI)
- Test: `tests/test_persistent_integration.py`

**Interfaces:**
- Consumes: `FatouGPWorker`, `WorkerDied` (Task 5).
- Produces:
  - `FatouGP` neue dataclass-Felder: `persistent: bool = True`, `init_timeout: float = 3600.0`, `eval_timeout: float = 600.0` (jeweils NACH den bestehenden Feldern, damit positionale Konstruktion nicht bricht).
  - `FatouGP.close() -> None`, `FatouGP.__enter__/__exit__`.
  - `FatouGP._spawn_worker(base: GPValue) -> FatouGPWorker` (Task 7 erweitert diese Methode um den Cache-Pfad).
  - Semantik: `_run_initialized` nutzt bei `persistent=True` einen gecachten Worker pro `_base_expr(base)`; bei `WorkerDied` genau ein Respawn-Retry; `persistent=False` = alter One-Shot-Pfad.

- [ ] **Step 1: Failing Tests schreiben**

`tests/test_persistent_integration.py`:

```python
from __future__ import annotations

import sys
import unittest
from pathlib import Path

from mpmath import mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend import FatouGP, find_default_fatou_gp, find_default_gp_exe


def make_gp(**kwargs) -> FatouGP:
    return FatouGP(gp_exe=find_default_gp_exe(), fatou_gp=find_default_fatou_gp(),
                   dps=50, nlim=20, nskip=4, looplim=30, **kwargs)


class PersistentIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        mp.dps = 80

    def test_persistent_matches_oneshot_real_base(self) -> None:
        with make_gp() as pers:
            oneshot = make_gp(persistent=False)
            xs = [mp.mpf("0.25"), mp.mpf("0.5")]
            a = pers.sexp_batch("e", xs)
            b = oneshot.sexp_batch("e", xs)
            for left, right in zip(a, b):
                self.assertLess(abs(left - right), mp.mpf("1e-30"))

    def test_persistent_matches_oneshot_complex_base(self) -> None:
        with make_gp() as pers:
            oneshot = make_gp(persistent=False)
            a = pers.sexp("1+I", mp.mpc("0.4", "0.2"))
            b = oneshot.sexp("1+I", mp.mpc("0.4", "0.2"))
            self.assertLess(abs(a - b), mp.mpf("1e-30"))

    def test_worker_is_reused_across_batches(self) -> None:
        with make_gp() as gp:
            gp.sexp("e", mp.mpf("0.25"))
            worker1 = gp._workers[gp._base_expr("e")]
            gp.slog("e", mp.mpf("2.0"))
            worker2 = gp._workers[gp._base_expr("e")]
            self.assertIs(worker1, worker2)

    def test_respawn_after_worker_death(self) -> None:
        with make_gp() as gp:
            gp.sexp("e", mp.mpf("0.25"))
            gp._workers[gp._base_expr("e")].close()  # simulate crash
            value = gp.sexp("e", mp.mpf("0.5"))
            self.assertLess(abs(value - mp.mpf("1.64635423375119458097")),
                            mp.mpf("1e-18"))

    def test_close_kills_all_workers(self) -> None:
        gp = make_gp()
        gp.sexp("e", mp.mpf("0.25"))
        workers = list(gp._workers.values())
        gp.close()
        self.assertTrue(all(not w.alive for w in workers))
        self.assertEqual(gp._workers, {})


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Tests laufen lassen, Fehlschlag verifizieren**

Run: `python -m unittest tests.test_persistent_integration -v`
Expected: ERROR (`unexpected keyword argument 'persistent'` bzw. `AttributeError: _workers`)

- [ ] **Step 3: Implementieren**

In `src/fatou_backend/gp_backend.py`:

Import ergänzen (nach `import mpmath as mp`):

```python
from .worker import FatouGPWorker, WorkerDied
```

`field` zum dataclass-Import ergänzen:

```python
from dataclasses import dataclass, field
```

`FatouGP`-Felder ergänzen (nach `quietmode: int = 1`):

```python
    persistent: bool = True
    init_timeout: float = 3600.0
    eval_timeout: float = 600.0
    _workers: dict = field(default_factory=dict, init=False, repr=False, compare=False)
```

Neue Methoden in `FatouGP` (nach `_base_expr`):

```python
    def _spawn_worker(self, base: GPValue) -> FatouGPWorker:
        return FatouGPWorker(
            self.gp_exe,
            self._init_lines(base),
            init_timeout=self.init_timeout,
            eval_timeout=self.eval_timeout,
        )

    def _worker_for(self, base: GPValue) -> FatouGPWorker:
        key = self._base_expr(base)
        worker = self._workers.get(key)
        if worker is None or not worker.alive:
            worker = self._spawn_worker(base)
            self._workers[key] = worker
        return worker

    def close(self) -> None:
        for worker in self._workers.values():
            worker.close()
        self._workers.clear()

    def __enter__(self) -> "FatouGP":
        return self

    def __exit__(self, *exc) -> None:
        self.close()
```

`_run_initialized` bekommt am Anfang das Routing (der Rest des bestehenden One-Shot-Codes bleibt unverändert darunter):

```python
    def _run_initialized(
        self,
        base: GPValue,
        expressions: list[str],
        digits: int | None = None,
    ) -> list[mp.mpc]:
        if digits is None:
            digits = max(self.dps - 8, 30)
        if self.persistent:
            try:
                return self._worker_for(base).eval(expressions)
            except WorkerDied:
                self._workers.pop(self._base_expr(base), None)
                return self._worker_for(base).eval(expressions)
        # legacy one-shot path: existing code from Task 1 continues unchanged here,
        # starting with:  lines = self._init_lines(base) + ['print("__BEGIN_RESULTS__")']
```

- [ ] **Step 4: Tests laufen lassen**

Run: `python -m unittest tests.test_persistent_integration -v`
Expected: PASS (5 Tests)

- [ ] **Step 5: Volle Suite (bestehende Tests laufen jetzt persistent!)**

Run: `python -m unittest discover -s tests -p "test_*.py" -v`
Expected: alle PASS — und deutlich schneller als vorher (die 12 Alt-Tests initialisieren Basen nur noch einmal pro `FatouGP`-Objekt). Sollte ein Alt-Test hängen/fehlschlagen, ist das ein Worker-Protokoll-Bug: mit `persistent=False` gegentesten, Ursache im Worker fixen — NICHT den Alt-Test ändern.

- [ ] **Step 6: M1-Messung**

Run: `python bench/benchmark.py --mode all --label m1-persistent`
Expected: `warm`-Zeilen jetzt um Größenordnungen schneller als `cold` (Init fällt weg); `min digits` unverändert ≥ 60 bei dps 80. Vergleich mit `bench/results/*-baseline-oneshot.json` dokumentieren:

Run: `python -c "import json,glob; b=json.load(open(sorted(glob.glob('bench/results/*baseline-oneshot.json'))[-1])); m=json.load(open(sorted(glob.glob('bench/results/*m1-persistent.json'))[-1])); bw={r['group']:r['digits_per_second'] for r in b['runs'] if r['mode']=='warm'}; mw={r['group']:r['digits_per_second'] for r in m['runs'] if r['mode']=='warm'}; [print(g, f'{mw[g]/bw[g]:.1f}x') for g in sorted(bw)]"`
Expected: Speedup-Faktoren pro Gruppe (Erwartung: >5× auf warm-Durchsatz; tatsächliche Zahl im Commit-Text festhalten).

- [ ] **Step 7: Commit**

```bash
git add src/fatou_backend/gp_backend.py tests/test_persistent_integration.py bench/results/
git commit -m "feat: persistente GP-Sessions als Default in FatouGP (M1); warm-Speedup <ZAHL>x ggue. One-Shot-Baseline

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

(`<ZAHL>` durch den gemessenen Median-Speedup aus Step 6 ersetzen.)

---

### Task 7: sexpinit-State-Cache (M2a)

**Files:**
- Create: `src/fatou_backend/state_cache.py`
- Modify: `src/fatou_backend/gp_backend.py` (`_spawn_worker` cache-aware, neues Feld `state_cache: bool = True`)
- Test: `tests/test_state_cache.py`

**Interfaces:**
- Consumes: `FatouGPWorker.eval_raw` (Task 5), `FatouGP._init_lines` (Task 1).
- Produces (`src/fatou_backend/state_cache.py`):
  - `cache_dir() -> Path` — `FATOU_CACHE_DIR`-Env oder `Path.home()/".cache"/"fatou_backend"`.
  - `cache_key(base_expr: str, dps: int, nlim: int, nskip: int, looplim: int, fatou_gp: Path, gp_exe: Path) -> str` — sha256-Hex über alle Felder inkl. Datei-Hash von fatou.gp und `gp_exe`-Pfadstring (32- vs 64-bit trennen).
  - `discover_state_var_names(worker: FatouGPWorker) -> list[str]` — via `\uv`.
  - `dump_state(worker, names, bin_path: Path) -> None` — `writebin` + Anker-Wert berechnen + JSON-Sidecar schreiben.
  - `restore_lines(names: list[str], bin_path: Path) -> list[str]` — GP-Zeilen für Restore.
  - `load_sidecar(bin_path: Path) -> dict | None` — `{"names": [...], "anchor_real": "...", "anchor_imag": "..."}` oder None wenn fehlt/korrupt.
  - `ANCHOR_EXPR = "sexp(0.5)"`.
- `FatouGP._spawn_worker` neue Semantik: Cache-Hit → Worker ohne `sexpinit` (read + restore), Anker validieren; Validierung scheitert → Cache-Dateien löschen, voller Spawn; Cache-Miss → voller Spawn, danach Discovery + Dump.

- [ ] **Step 1: Failing Tests schreiben**

`tests/test_state_cache.py`:

```python
from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from mpmath import mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend import FatouGP, find_default_fatou_gp, find_default_gp_exe
from fatou_backend import state_cache


def make_gp(**kwargs) -> FatouGP:
    return FatouGP(gp_exe=find_default_gp_exe(), fatou_gp=find_default_fatou_gp(),
                   dps=50, nlim=20, nskip=4, looplim=30, **kwargs)


class StateCacheTests(unittest.TestCase):
    def setUp(self) -> None:
        mp.dps = 80
        self._tmp = tempfile.TemporaryDirectory()
        self._env = mock.patch.dict(os.environ, {"FATOU_CACHE_DIR": self._tmp.name})
        self._env.start()

    def tearDown(self) -> None:
        self._env.stop()
        self._tmp.cleanup()

    def test_cache_key_depends_on_dps(self) -> None:
        gp = make_gp()
        k1 = state_cache.cache_key("exp(1)", 50, 20, 4, 30, gp.fatou_gp, gp.gp_exe)
        k2 = state_cache.cache_key("exp(1)", 60, 20, 4, 30, gp.fatou_gp, gp.gp_exe)
        self.assertNotEqual(k1, k2)
        self.assertEqual(len(k1), 64)

    def test_cache_written_on_first_use_and_used_on_second(self) -> None:
        with make_gp() as gp1:
            v1 = gp1.sexp("e", mp.mpf("0.5"))
        bins = list(Path(self._tmp.name).glob("*.gpbin"))
        self.assertEqual(len(bins), 1)
        self.assertTrue(bins[0].with_suffix(".json").exists())
        with make_gp() as gp2:
            v2 = gp2.sexp("e", mp.mpf("0.5"))
            # cached worker must not have run sexpinit
            worker = gp2._workers[gp2._base_expr("e")]
            self.assertFalse(any("sexpinit" in ln for ln in worker.debug_init_lines))
        self.assertLess(abs(v1 - v2), mp.mpf("1e-30"))

    def test_corrupt_cache_falls_back_to_full_init(self) -> None:
        with make_gp() as gp1:
            v1 = gp1.sexp("e", mp.mpf("0.5"))
        bin_path = next(Path(self._tmp.name).glob("*.gpbin"))
        bin_path.write_bytes(b"garbage")
        with make_gp() as gp2:
            v2 = gp2.sexp("e", mp.mpf("0.5"))
        self.assertLess(abs(v1 - v2), mp.mpf("1e-30"))

    def test_state_cache_disabled(self) -> None:
        with make_gp(state_cache=False) as gp:
            gp.sexp("e", mp.mpf("0.5"))
        self.assertEqual(list(Path(self._tmp.name).glob("*")), [])


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Tests laufen lassen, Fehlschlag verifizieren**

Run: `python -m unittest tests.test_state_cache -v`
Expected: ERROR `ImportError: cannot import name 'state_cache'`

- [ ] **Step 3: `state_cache.py` implementieren**

```python
"""Persist sexpinit state per (base, dps, knobs, fatou-hash) via writebin/read.

Discovery uses the \\uv metacommand (verified on gp 2.17.3): it lists all user
variables as "name =" lines with indented value lines. We dump every variable
into one binary vector; functions are NOT dumped (re-reading fatou.gp is cheap
and restores them).
"""
from __future__ import annotations

import hashlib
import json
import os
import re
from pathlib import Path

import mpmath as mp

ANCHOR_EXPR = "sexp(0.5)"
_NAME_RE = re.compile(r"^([A-Za-z_][A-Za-z0-9_]*) =$")
_EXCLUDED = {"vv"}  # our own eval temp variable


def cache_dir() -> Path:
    env = os.getenv("FATOU_CACHE_DIR")
    base = Path(env) if env else Path.home() / ".cache" / "fatou_backend"
    base.mkdir(parents=True, exist_ok=True)
    return base


def cache_key(base_expr: str, dps: int, nlim: int, nskip: int, looplim: int,
              fatou_gp: Path, gp_exe: Path) -> str:
    fatou_hash = hashlib.sha256(Path(fatou_gp).read_bytes()).hexdigest()
    blob = "|".join([base_expr, str(dps), str(nlim), str(nskip), str(looplim),
                     fatou_hash, str(gp_exe)])
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def discover_state_var_names(worker) -> list[str]:
    output = worker.eval_raw(["\\uv"], "__UVDONE__", timeout=120.0)
    names = []
    for line in output:
        match = _NAME_RE.match(line)
        if match and match.group(1) not in _EXCLUDED:
            names.append(match.group(1))
    return names


def dump_state(worker, names: list[str], bin_path: Path) -> None:
    gp_path = bin_path.as_posix()
    worker.eval_raw([f'writebin("{gp_path}", [{", ".join(names)}]);'],
                    "__DUMPDONE__", timeout=300.0)
    anchor = worker.eval([ANCHOR_EXPR])[0]
    sidecar = {
        "names": names,
        "anchor_real": mp.nstr(mp.re(anchor), 40, min_fixed=0, max_fixed=0),
        "anchor_imag": mp.nstr(mp.im(anchor), 40, min_fixed=0, max_fixed=0),
    }
    bin_path.with_suffix(".json").write_text(json.dumps(sidecar), encoding="utf-8")


def restore_lines(names: list[str], bin_path: Path) -> list[str]:
    gp_path = bin_path.as_posix()
    lines = [f'__fst = read("{gp_path}");']
    for idx, name in enumerate(names, start=1):
        lines.append(f"{name} = __fst[{idx}];")
    return lines


def load_sidecar(bin_path: Path) -> dict | None:
    sidecar_path = bin_path.with_suffix(".json")
    if not bin_path.exists() or not sidecar_path.exists():
        return None
    try:
        data = json.loads(sidecar_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(data.get("names"), list) or "anchor_real" not in data:
        return None
    return data


def drop_cache(bin_path: Path) -> None:
    bin_path.unlink(missing_ok=True)
    bin_path.with_suffix(".json").unlink(missing_ok=True)
```

- [ ] **Step 4: `FatouGP._spawn_worker` cache-aware machen**

In `src/fatou_backend/gp_backend.py`: Import `from . import state_cache` ergänzen, Feld `state_cache: bool = True` zur dataclass (nach `eval_timeout`). `FatouGPWorker` braucht das Attribut `debug_init_lines` (die Init-Zeilen, mit denen er gestartet wurde) — in `worker.py` `__init__` als erste Zeile `self.debug_init_lines = list(init_lines)` ergänzen.

`_spawn_worker` ersetzen:

```python
    def _spawn_worker(self, base: GPValue) -> FatouGPWorker:
        if not self.state_cache:
            return self._spawn_full(base)
        key = state_cache.cache_key(
            self._base_expr(base), self.dps, self.nlim, self.nskip,
            self.looplim, self.fatou_gp, self.gp_exe)
        bin_path = state_cache.cache_dir() / f"{key}.gpbin"
        sidecar = state_cache.load_sidecar(bin_path)
        if sidecar is not None:
            worker = self._try_spawn_cached(base, bin_path, sidecar)
            if worker is not None:
                return worker
            state_cache.drop_cache(bin_path)
        worker = self._spawn_full(base)
        try:
            names = state_cache.discover_state_var_names(worker)
            state_cache.dump_state(worker, names, bin_path)
        except (RuntimeError, OSError):
            state_cache.drop_cache(bin_path)  # cache is best-effort
        return worker

    def _spawn_full(self, base: GPValue) -> FatouGPWorker:
        return FatouGPWorker(
            self.gp_exe, self._init_lines(base),
            init_timeout=self.init_timeout, eval_timeout=self.eval_timeout)

    def _try_spawn_cached(self, base: GPValue, bin_path, sidecar) -> FatouGPWorker | None:
        init = self._init_lines(base)
        lines = init[:3] + state_cache.restore_lines(sidecar["names"], bin_path)
        try:
            worker = FatouGPWorker(
                self.gp_exe, lines,
                init_timeout=self.init_timeout, eval_timeout=self.eval_timeout)
            anchor = worker.eval([state_cache.ANCHOR_EXPR])[0]
        except (RuntimeError, OSError):
            return None
        expected = mp.mpc(mp.mpf(sidecar["anchor_real"]), mp.mpf(sidecar["anchor_imag"]))
        tolerance = mp.mpf(10) ** (-(min(self.dps, 35)))
        if abs(anchor - expected) > tolerance * max(1, abs(expected)):
            worker.close()
            return None
        return worker
```

Anmerkung: `init[:3]` = precision, `read(fatou.gp)`, `quietmode` — bewusst OHNE die `sexpinit`-Zeile (Index 3). Der Sidecar-Anker ist mit 40 Stellen gespeichert; die Toleranz `10^-min(dps,35)` bleibt deshalb unterhalb der Sidecar-Genauigkeit und ist streng genug, um jeden echten State-Schaden zu erkennen.

- [ ] **Step 5: Tests laufen lassen**

Run: `python -m unittest tests.test_state_cache -v`
Expected: PASS (4 Tests)

- [ ] **Step 6: Volle Suite + Cache-Verhalten im Benchmark messen**

Run: `python -m unittest discover -s tests -p "test_*.py" -v`
Expected: alle PASS. Achtung: die Alt-Tests schreiben jetzt in den ECHTEN Cache (`~/.cache/fatou_backend`) — das ist gewollt (zweiter Suite-Lauf wird schneller). Prüfen: `python -m unittest discover -s tests -p "test_gp_backend.py"` zweimal laufen lassen; zweiter Lauf messbar schneller.

Run: `python bench/benchmark.py --mode all --label m2-cache`
Expected: `cold`-Zeiten beim ZWEITEN Lauf des Kommandos deutlich unter dem M1-Stand (cold enthält jetzt nur noch Prozessstart + read + restore statt sexpinit). Beide Läufe ausführen, zweiten als Ergebnis nehmen.

- [ ] **Step 7: Commit**

```bash
git add src/fatou_backend/state_cache.py src/fatou_backend/gp_backend.py src/fatou_backend/worker.py tests/test_state_cache.py bench/results/
git commit -m "feat: sexpinit-State-Cache via writebin/read mit Anker-Validierung (M2a)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 8: FatouGPPool + mixed_phase-Anbindung (M2b)

**Files:**
- Create: `src/fatou_backend/pool.py`
- Modify: `src/fatou_backend/gp_backend.py` (Feld `n_workers: int = 1`, Routing in `_run_initialized`, `close()` schließt Pools)
- Modify: `src/fatou_backend/mixed_phase.py` (`--workers`-Flag)
- Test: `tests/test_pool.py`

**Interfaces:**
- Consumes: `FatouGP._spawn_worker` (Task 7), `FatouGPWorker.eval`.
- Produces:
  - `class FatouGPPool` in `pool.py`:
    - `__init__(self, spawn: Callable[[], FatouGPWorker], n_workers: int)`
    - `eval(self, expressions: list[str]) -> list[mp.mpc]` — Round-Robin-Chunks (`expressions[i::n]`), parallel via `ThreadPoolExecutor`, Ergebnisse in Originalreihenfolge.
    - `close(self) -> None`
  - `FatouGP(n_workers=4)` routet Batches mit ≥ 2·n_workers Ausdrücken über den Pool (kleinere Batches: einzelner Worker, Pool-Overhead lohnt nicht).
  - `mixed_phase --workers N` → `FatouGP(n_workers=N)`.

- [ ] **Step 1: Failing Tests schreiben**

`tests/test_pool.py`:

```python
from __future__ import annotations

import sys
import unittest
from pathlib import Path

from mpmath import mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend import FatouGP, find_default_fatou_gp, find_default_gp_exe
from fatou_backend.pool import FatouGPPool


def make_gp(**kwargs) -> FatouGP:
    return FatouGP(gp_exe=find_default_gp_exe(), fatou_gp=find_default_fatou_gp(),
                   dps=50, nlim=20, nskip=4, looplim=30, **kwargs)


class PoolTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        mp.dps = 80

    def test_pool_preserves_order_and_values(self) -> None:
        with make_gp() as single:
            xs = [mp.mpf(j) / 8 for j in range(8)]
            want = single.sexp_batch("e", xs)
        with make_gp(n_workers=2) as pooled:
            got = pooled.sexp_batch("e", xs)
        self.assertEqual(len(got), len(want))
        for left, right in zip(got, want):
            self.assertLess(abs(left - right), mp.mpf("1e-30"))

    def test_pool_class_round_robin(self) -> None:
        gp = make_gp()
        pool = FatouGPPool(lambda: gp._spawn_worker("e"), n_workers=2)
        try:
            values = pool.eval([f"sexp({j}/8)" for j in range(5)])
            self.assertEqual(len(values), 5)
            direct = gp._spawn_worker("e")
            try:
                expected = direct.eval([f"sexp({j}/8)" for j in range(5)])
            finally:
                direct.close()
            for left, right in zip(values, expected):
                self.assertLess(abs(left - right), mp.mpf("1e-30"))
        finally:
            pool.close()
            gp.close()

    def test_small_batch_avoids_pool(self) -> None:
        with make_gp(n_workers=4) as gp:
            gp.sexp("e", mp.mpf("0.5"))   # 1 expression < 2*4 -> single worker
            self.assertEqual(gp._pools, {})


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Tests laufen lassen, Fehlschlag verifizieren**

Run: `python -m unittest tests.test_pool -v`
Expected: ERROR `ModuleNotFoundError: No module named 'fatou_backend.pool'`

- [ ] **Step 3: Implementieren**

`src/fatou_backend/pool.py`:

```python
"""Round-robin pool of persistent workers for large batches."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from typing import Callable

import mpmath as mp

from .worker import FatouGPWorker


class FatouGPPool:
    def __init__(self, spawn: Callable[[], FatouGPWorker], n_workers: int) -> None:
        if n_workers < 1:
            raise ValueError("n_workers must be >= 1")
        self._workers = [spawn() for _ in range(n_workers)]
        self._executor = ThreadPoolExecutor(max_workers=n_workers)

    def eval(self, expressions: list[str]) -> list[mp.mpc]:
        n = len(self._workers)
        chunks = [expressions[i::n] for i in range(n)]
        futures = [
            self._executor.submit(worker.eval, chunk)
            for worker, chunk in zip(self._workers, chunks) if chunk
        ]
        chunk_results = [f.result() for f in futures]
        merged: list[mp.mpc | None] = [None] * len(expressions)
        for i, results in enumerate(chunk_results):
            for j, value in enumerate(results):
                merged[i + j * n] = value
        return merged  # type: ignore[return-value]

    def close(self) -> None:
        for worker in self._workers:
            worker.close()
        self._executor.shutdown(wait=False)
```

In `gp_backend.py`: Import `from .pool import FatouGPPool`; Felder ergänzen:

```python
    n_workers: int = 1
    _pools: dict = field(default_factory=dict, init=False, repr=False, compare=False)
```

Routing in `_run_initialized` (den `if self.persistent:`-Block ersetzen):

```python
        if self.persistent:
            if self.n_workers > 1 and len(expressions) >= 2 * self.n_workers:
                return self._pool_for(base).eval(expressions)
            try:
                return self._worker_for(base).eval(expressions)
            except WorkerDied:
                self._workers.pop(self._base_expr(base), None)
                return self._worker_for(base).eval(expressions)
```

Neue Methode + `close()` erweitern:

```python
    def _pool_for(self, base: GPValue) -> FatouGPPool:
        key = self._base_expr(base)
        pool = self._pools.get(key)
        if pool is None:
            pool = FatouGPPool(lambda: self._spawn_worker(base), self.n_workers)
            self._pools[key] = pool
        return pool

    def close(self) -> None:
        for worker in self._workers.values():
            worker.close()
        self._workers.clear()
        for pool in self._pools.values():
            pool.close()
        self._pools.clear()
```

In `mixed_phase.py`: argparse-Zeile ergänzen (`parse_args`):

```python
    parser.add_argument("--workers", type=int, default=1)
```

und in `main()` den FatouGP-Aufruf erweitern:

```python
    gp = FatouGP(gp_exe=args.gp_exe, fatou_gp=args.fatou_gp, dps=args.dps,
                 looplim=max(35, args.dps - 20), n_workers=args.workers)
```

- [ ] **Step 4: Tests laufen lassen**

Run: `python -m unittest tests.test_pool -v`
Expected: PASS (3 Tests)

- [ ] **Step 5: Volle Suite + M2-Messung**

Run: `python -m unittest discover -s tests -p "test_*.py" -v`
Expected: alle PASS

M2-Messung (PowerShell, Wandzeit vergleichen):

```powershell
Measure-Command { python -m fatou_backend.mixed_phase --outer e --inner 2 --samples 64 --workers 1 } | Select-Object TotalSeconds
Measure-Command { python -m fatou_backend.mixed_phase --outer e --inner 2 --samples 64 --workers 4 } | Select-Object TotalSeconds
```

Expected: `--workers 4` spürbar schneller auf den sexp/slog-Batches (Faktor ~2–3 bei 4 Workern; der sexpinit-Cache aus Task 7 macht die zusätzlichen Worker-Starts billig). Beide Zeiten im Commit-Text notieren.

- [ ] **Step 6: Commit**

```bash
git add src/fatou_backend/pool.py src/fatou_backend/gp_backend.py src/fatou_backend/mixed_phase.py tests/test_pool.py
git commit -m "feat: FatouGPPool (Round-Robin ueber persistente Worker) + mixed_phase --workers (M2b)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 9: research/-Struktur — Fork, Gate, program.md, journal.md (M3a)

**Files:**
- Create: `src/fatou_backend/vendor/fatou_fork.gp` (byte-identische Kopie)
- Create: `research/gate.py`, `research/program.md`, `research/journal.md`
- Test: `tests/test_gate.py`

**Interfaces:**
- Consumes: `bench.workload`, `bench.metrics`, `research/reference/values.json`, `FatouGP`.
- Produces:
  - `research/gate.py` CLI: `python research/gate.py [--fork PATH]` (Default-Fork: `src/fatou_backend/vendor/fatou_fork.gp`). Exit 0 = PASS, Exit 1 = FAIL; druckt pro Gruppe `PASS/FAIL`-Zeilen.
  - Funktion `run_gate(fork_path: Path) -> tuple[bool, list[str]]` (für Tests importierbar).
  - Gate-Kriterien (eingefroren): (a) jeder Workload-Case stimmt mit der Referenz auf ≥ `case.dps - 5` Stellen überein (`correct_digits`-Metrik, cap = dps), (b) Roundtrip `sexp(slog(y)) - y` für `y ∈ {1.5, 2.5, 3.5}` je Basis e/2/10 bei dps 80 hat `|residual| < 1e-40`, (c) komplexe Anker-Cases aus dem Workload bestehen (a) mit.
- **Ab diesem Task sind `research/gate.py` und `research/reference/` EINGEFROREN.**

- [ ] **Step 1: Fork anlegen**

```bash
cd "C:\Users\janis\Documents\Tetration\gp-tetration"
cp src/fatou_backend/vendor/fatou.gp src/fatou_backend/vendor/fatou_fork.gp
git add src/fatou_backend/vendor/fatou_fork.gp
```

- [ ] **Step 2: Failing Tests schreiben**

`tests/test_gate.py`:

```python
from __future__ import annotations

import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "research"))

from gate import run_gate

REPO = Path(__file__).resolve().parents[1]
ORIGINAL = REPO / "src" / "fatou_backend" / "vendor" / "fatou.gp"
FORK = REPO / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"

RUN_SLOW = os.getenv("FATOU_BACKEND_RUN_SLOW") == "1"


class GateTests(unittest.TestCase):
    def test_fork_starts_byte_identical(self) -> None:
        self.assertEqual(ORIGINAL.read_bytes(), FORK.read_bytes())

    @unittest.skipUnless(RUN_SLOW, "Set FATOU_BACKEND_RUN_SLOW=1")
    def test_gate_passes_on_unmodified_fork(self) -> None:
        ok, report = run_gate(FORK)
        self.assertTrue(ok, "\n".join(report))

    @unittest.skipUnless(RUN_SLOW, "Set FATOU_BACKEND_RUN_SLOW=1")
    def test_gate_fails_on_sabotaged_fork(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sabotaged = Path(tmp) / "fatou_sabotaged.gp"
            shutil.copy(FORK, sabotaged)
            with sabotaged.open("a", encoding="utf-8") as handle:
                handle.write("\nslog(z) = { abel(z*lnb+k-1)+rslog + 1e-20; }\n")
            ok, report = run_gate(sabotaged)
            self.assertFalse(ok)
            self.assertTrue(any("FAIL" in line for line in report))


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 3: Tests laufen lassen, Fehlschlag verifizieren**

Run: `python -m unittest tests.test_gate -v`
Expected: ERROR `ModuleNotFoundError: No module named 'gate'`

- [ ] **Step 4: `research/gate.py` implementieren**

```python
"""Immutable correctness gate for fatou_fork.gp experiments.

FROZEN after task 9 of the 2026-07-09 plan. Loop experiments must never edit
this file or research/reference/. If the gate itself needs a change, that is a
human decision outside the loop, documented in the journal.
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import mpmath as mp

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

from bench.metrics import correct_digits, load_reference
from bench.workload import all_cases, case_id
from fatou_backend import FatouGP, find_default_gp_exe

REFERENCE_PATH = REPO / "research" / "reference" / "values.json"
DEFAULT_FORK = REPO / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
AGREE_MARGIN = 5          # required digits = case.dps - AGREE_MARGIN
ROUNDTRIP_MAX = mp.mpf("1e-40")
ROUNDTRIP_BASES = ("e", "2", "10")
ROUNDTRIP_YS = ("1.5", "2.5", "3.5")
ROUNDTRIP_DPS = 80


def run_gate(fork_path: Path) -> tuple[bool, list[str]]:
    gp_exe = find_default_gp_exe()
    report: list[str] = []
    ok = True

    groups: dict[tuple[str, int], list] = defaultdict(list)
    for c in all_cases(deep=False):
        groups[(c.base, c.dps)].append(c)

    for (base, dps), cases in sorted(groups.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        mp.mp.dps = dps + 40
        reference = load_reference(REFERENCE_PATH)
        with FatouGP(gp_exe=gp_exe, fatou_gp=fork_path, dps=dps,
                     looplim=max(35, dps - 20), state_cache=False) as gp:
            try:
                results = gp.eval_batch(base, [f"{c.kind}({c.arg})" for c in cases])
            except RuntimeError as exc:
                ok = False
                report.append(f"FAIL {base}|{dps}: gp error: {exc}")
                continue
        worst = min(correct_digits(v, reference[case_id(c)], cap=dps)
                    for c, v in zip(cases, results))
        required = dps - AGREE_MARGIN
        line_ok = worst >= required
        ok = ok and line_ok
        report.append(f"{'PASS' if line_ok else 'FAIL'} agree {base}|{dps}: "
                      f"worst {worst:.1f} digits (required {required})")

    mp.mp.dps = ROUNDTRIP_DPS + 40
    for base in ROUNDTRIP_BASES:
        with FatouGP(gp_exe=gp_exe, fatou_gp=fork_path, dps=ROUNDTRIP_DPS,
                     state_cache=False) as gp:
            try:
                residuals = gp.roundtrip_residuals(base, [mp.mpf(y) for y in ROUNDTRIP_YS])
            except RuntimeError as exc:
                ok = False
                report.append(f"FAIL roundtrip {base}: gp error: {exc}")
                continue
        worst_res = max(abs(r) for r in residuals)
        line_ok = worst_res < ROUNDTRIP_MAX
        ok = ok and line_ok
        report.append(f"{'PASS' if line_ok else 'FAIL'} roundtrip {base}: "
                      f"worst {mp.nstr(worst_res, 6)} (max {mp.nstr(ROUNDTRIP_MAX, 3)})")

    return ok, report


def main() -> None:
    parser = argparse.ArgumentParser(description="fatou_fork.gp correctness gate")
    parser.add_argument("--fork", default=str(DEFAULT_FORK))
    args = parser.parse_args()
    ok, report = run_gate(Path(args.fork))
    print("\n".join(report))
    print("GATE:", "PASS" if ok else "FAIL")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
```

Anmerkung: `state_cache=False` im Gate — der Cache ist nach `fatou_gp`-Datei-Hash geschlüsselt und wäre korrekt, aber das Gate soll den Fork kompromisslos frisch initialisieren (keine Cache-Pfad-Interaktion im Urteil).

- [ ] **Step 5: Tests laufen lassen (Gate-Läufe sind slow-guarded)**

Run: `FATOU_BACKEND_RUN_SLOW=1 python -m unittest tests.test_gate -v` (Bash; PowerShell: `$env:FATOU_BACKEND_RUN_SLOW="1"; python -m unittest tests.test_gate -v`)
Expected: PASS (3 Tests). Laufzeit: das Gate initialisiert 8 (basis, dps)-Gruppen + 3 Roundtrip-Basen, und das zweimal (pass + sabotage) — zweistellige Minutenzahl ist ok; deshalb der `FATOU_BACKEND_RUN_SLOW`-Guard (Repo-Konvention aus `test_gp_backend_slow.py`), damit die Schnell-Suite unter 2 Minuten bleibt. Falls `test_gate_passes_on_unmodified_fork` an der Roundtrip-Schwelle scheitert: tatsächlichen worst-Wert aus dem Report ablesen; liegt er über `1e-40`, ist `ROUNDTRIP_MAX` in gate.py EINMALIG vor dem Einfrieren auf `worst * 100` anzuheben (im Commit dokumentieren). Danach gilt die Datei als eingefroren.

- [ ] **Step 6: `research/program.md` schreiben**

```markdown
# Autoresearch-Programm: fatou_fork.gp beschleunigen

## Ziel

Verifizierte korrekte Stellen pro Sekunde maximieren (bench/benchmark.py auf dem
eingefrorenen Workload), ohne das Gate zu verletzen.

## Unantastbar (Anti-Gaming)

- `research/gate.py`, `research/reference/values.json`, `bench/workload.py`,
  `bench/metrics.py`: NIEMALS aendern. Wenn eine Aenderung dort noetig scheint,
  ist das ein Mensch-Entscheid ausserhalb des Loops -> Journal-Eintrag, stoppen.
- Referenz ist immer die Original-`fatou.gp`; sie wird nie editiert.

## Experiment-Protokoll (ein Experiment = eine Mutation)

1. Hypothese im Journal notieren (was wird geaendert, warum sollte es schneller sein).
2. Genau EINE Mutation an `src/fatou_backend/vendor/fatou_fork.gp` ODER an den
   Wrapper-Knobs (nlim/nskip/looplim/dps-Staffelung) vornehmen.
3. `python research/gate.py` -> bei FAIL: revert, Journal-Eintrag, naechstes Experiment.
4. Bei PASS: `python bench/benchmark.py --fatou src/fatou_backend/vendor/fatou_fork.gp
   --label exp-<NNN> --repeat 3` (warm-Median zaehlt).
5. Behalten nur, wenn die Metrik um mehr als die Rausch-Schwelle aus
   `research/noise.json` steigt; sonst revert.
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
```

- [ ] **Step 7: `research/journal.md` anlegen**

```markdown
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

## exp-000 (2026-07-09) — Baseline
- **Hypothese:** — (Ausgangszustand nach M1/M2)
- **Mutation:** keine; fatou_fork.gp == fatou.gp
- **Gate:** PASS (siehe Task-9-Testlauf)
- **Benchmark:** siehe bench/results/ (Labels baseline-oneshot, m1-persistent, m2-cache)
- **Entscheidung:** Referenzpunkt
- **Learnings:** Speedup-Treppe One-Shot -> persistent -> Cache dokumentiert.
```

- [ ] **Step 8: Volle Suite laufen lassen**

Run: `python -m unittest discover -s tests -p "test_*.py" -v`
Expected: alle PASS

- [ ] **Step 9: Commit**

```bash
git add src/fatou_backend/vendor/fatou_fork.gp research/gate.py research/program.md research/journal.md tests/test_gate.py
git commit -m "feat: research-Struktur - fatou_fork.gp, eingefrorenes Gate, program.md, Journal (M3a)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

### Task 10: Rausch-Kalibrierung + README (M3b)

**Files:**
- Create: `research/calibrate_noise.py`, `research/noise.json` (generiert)
- Modify: `README.md` (Abschnitte zu persistenten Sessions, Cache, Pool, Benchmark, research-Loop)
- Test: `tests/test_calibrate_noise.py`

**Interfaces:**
- Consumes: Benchmark-Ergebnis-JSON-Schema (Task 4).
- Produces:
  - `research.calibrate_noise.noise_threshold_pct(samples: list[float]) -> float` — `max(3.0, 2 * 100 * (max-min)/median)`.
  - `research/noise.json`: `{"threshold_pct": <float>, "samples": [...], "source": "<results-json>"}`.

- [ ] **Step 1: Failing Test schreiben**

`tests/test_calibrate_noise.py`:

```python
from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "research"))

from calibrate_noise import noise_threshold_pct


class NoiseThresholdTests(unittest.TestCase):
    def test_floor_is_three_percent(self) -> None:
        self.assertEqual(noise_threshold_pct([1000.0, 1000.1, 999.9]), 3.0)

    def test_wide_spread_raises_threshold(self) -> None:
        # spread (1200-800)/1000 = 40% -> threshold 80%
        self.assertAlmostEqual(noise_threshold_pct([800.0, 1000.0, 1200.0]), 80.0)

    def test_requires_at_least_three_samples(self) -> None:
        with self.assertRaises(ValueError):
            noise_threshold_pct([1.0, 2.0])


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Test laufen lassen, Fehlschlag verifizieren**

Run: `python -m unittest tests.test_calibrate_noise -v`
Expected: ERROR `ModuleNotFoundError: No module named 'calibrate_noise'`

- [ ] **Step 3: Implementieren**

`research/calibrate_noise.py`:

```python
"""Calibrate the benchmark noise threshold from repeated warm runs."""
from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
NOISE_PATH = REPO / "research" / "noise.json"


def noise_threshold_pct(samples: list[float]) -> float:
    if len(samples) < 3:
        raise ValueError("need at least 3 samples")
    median = statistics.median(samples)
    spread_pct = 100.0 * (max(samples) - min(samples)) / median
    return max(3.0, 2.0 * spread_pct)


def main() -> None:
    parser = argparse.ArgumentParser(description="calibrate benchmark noise")
    parser.add_argument("--repeat", type=int, default=5)
    args = parser.parse_args()

    label = "noise-calibration"
    subprocess.run(
        [sys.executable, str(REPO / "bench" / "benchmark.py"),
         "--mode", "throughput", "--repeat", str(args.repeat), "--label", label],
        check=True,
    )
    results = sorted((REPO / "bench" / "results").glob(f"*-{label}.json"))
    payload = json.loads(results[-1].read_text(encoding="utf-8"))
    # total warm digits/s per repetition, summed over groups
    per_rep: dict[int, float] = {}
    for run in payload["runs"]:
        if run["mode"] == "warm":
            per_rep[run["rep"]] = per_rep.get(run["rep"], 0.0) + run["digits_per_second"]
    samples = [per_rep[k] for k in sorted(per_rep)]
    threshold = noise_threshold_pct(samples)
    NOISE_PATH.write_text(
        json.dumps({"threshold_pct": threshold, "samples": samples,
                    "source": results[-1].name}, indent=1),
        encoding="utf-8")
    print(f"[calibrate_noise] threshold {threshold:.1f}% from {len(samples)} runs -> {NOISE_PATH}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Tests laufen lassen**

Run: `python -m unittest tests.test_calibrate_noise -v`
Expected: PASS (3 Tests)

- [ ] **Step 5: Kalibrierung ausführen**

Run: `python research/calibrate_noise.py --repeat 5`
Expected: Benchmark läuft 5× (warm dominiert die Kosten dank Cache), `research/noise.json` mit `threshold_pct` ≥ 3.0. Bei `threshold_pct` > 15 % ist die Maschine zu unruhig — Hintergrundlast schließen und wiederholen, bevor Experimente starten.

- [ ] **Step 6: README ergänzen**

In `README.md` nach dem Abschnitt „## Tests" anfügen:

```markdown
## Performance-Architektur

- Persistente GP-Sessions sind Default (`FatouGP(persistent=False)` fuer den
  alten One-Shot-Pfad). `gp.close()` oder Context-Manager beendet die Worker.
- Der sexpinit-State wird pro (Basis, dps, Knobs, fatou-Hash) in
  `~/.cache/fatou_backend/` gecacht (`FATOU_CACHE_DIR` uebersteuert;
  `FatouGP(state_cache=False)` schaltet ab).
- `FatouGP(n_workers=N)` bzw. `mixed_phase --workers N` parallelisiert grosse
  Batches ueber N Worker-Prozesse.

## Benchmark & Autoresearch-Loop

- `python bench/benchmark.py --mode all --label <name>` misst verifizierte
  korrekte Stellen pro Sekunde gegen `research/reference/values.json`.
- `python research/gate.py` prueft `fatou_fork.gp` gegen die eingefrorene
  Referenz (Exit 0 = PASS).
- Loop-Regeln: `research/program.md`; Historie: `research/journal.md`.
```

- [ ] **Step 7: Volle Suite + Schlussmessung**

Run: `python -m unittest discover -s tests -p "test_*.py" -v`
Expected: alle PASS

Run: `python bench/benchmark.py --mode all --repeat 3 --label m3-final`
Expected: konsistente warm-Werte (Spread unter der kalibrierten Schwelle); Ergebnis als Schlussstand von M1–M3 in `bench/results/`.

- [ ] **Step 8: Commit**

```bash
git add research/calibrate_noise.py research/noise.json tests/test_calibrate_noise.py README.md bench/results/
git commit -m "feat: Rausch-Kalibrierung + README-Doku; M3 abgeschlossen, Loop startklar

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

---

## Nacharbeiten (außerhalb dieses Plans)

- **Deep-Referenzen (dps 500/1000):** nach M2-Speedups `python bench/make_reference.py --deep` laufen lassen (ergänzt `values.json`; das ist die EINZIGE erlaubte nachträgliche Änderung an `research/reference/` — sie fügt Werte hinzu, ändert keine bestehenden), committen, dann in Benchmark-Läufen `--deep` nutzen.
- **Loop-Betrieb (M4):** Experimente nach `research/program.md` fahren — manuell oder autonom via `/loop`. Erster Kandidat laut fatou.gp-Header: `limitp=16` („default precision can be set to 16 to speed up calculations").
