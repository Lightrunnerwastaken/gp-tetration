"""A/B cold-init timing for engine mutations.

Times sexpinit() for one or more .gp engine files at one or more dps, and
prints the resulting anchor sexp(0.5) so accuracy and speed come from the
same run. No state cache is involved (plain GP), so every number is a cold
initialization -- which is the cost block the exponent lives in.

    python research/tools/ab_init.py --dps 200 300 \
        --engine src/fatou_backend/vendor/fatou_fork.gp \
        --engine research/tools/_exp041.gp --repeat 2

Agreement between two engines is reported as the number of matching leading
decimal digits of the anchor.
"""
from __future__ import annotations

import argparse
import re as _re
import subprocess
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DEFAULT_GP = r"C:\Program Files (x86)\Pari64-2-17-3\gp.exe"


def _script(engine: Path, base: str, dps: int) -> str:
    return "\n".join([
        "default(parisizemax, 8589934592);",
        f"default(realprecision, {dps});",
        f'read("{engine.as_posix()}");',
        f"default(realprecision, {dps});",
        "quietmode=1;",
        f"abbase = {base};",
        "gettime();",
        "sexpinit(abbase,0,0,0);",
        'print("AB-MS ", gettime());',
        'print("AB-VAL ", sexp(0.5));',
        "quit;",
    ]) + "\n"


def run(engine: Path, base: str, dps: int, gp: str) -> tuple[int, str, float]:
    t0 = time.time()
    proc = subprocess.run(
        [gp, "-q", "-f"], input=_script(engine, base, dps),
        capture_output=True, text=True, encoding="utf-8", errors="replace",
        cwd=str(REPO),
    )
    wall = time.time() - t0
    ms = _re.search(r"AB-MS\s+(\d+)", proc.stdout)
    val = _re.search(r"AB-VAL\s+([-\d.eE]+)", proc.stdout)
    if not ms or not val:
        raise SystemExit(f"no result for {engine.name} @dps {dps}\n"
                         f"{proc.stdout[-1500:]}\n{proc.stderr[-800:]}")
    return int(ms.group(1)), val.group(1), wall


def agree_digits(a: str, b: str) -> int:
    n = 0
    for ca, cb in zip(a, b):
        if ca != cb:
            break
        if ca.isdigit():
            n += 1
    return n


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", action="append", required=True)
    ap.add_argument("--dps", type=int, nargs="+", default=[200])
    ap.add_argument("--base", default="exp(1)")
    ap.add_argument("--repeat", type=int, default=1)
    ap.add_argument("--gp", default=DEFAULT_GP)
    args = ap.parse_args()

    engines = [Path(e) if Path(e).is_absolute() else REPO / e for e in args.engine]
    for dps in args.dps:
        results = []
        for eng in engines:
            best, val = None, None
            for _ in range(args.repeat):
                ms, v, wall = run(eng, args.base, dps, args.gp)
                best = ms if best is None else min(best, ms)
                val = v
            results.append((eng.name, best, val))
            print(f"  dps {dps:5d}  {eng.name:28s} {best/1000:9.2f}s")
        if len(results) == 2:
            (_, m0, v0), (_, m1, v1) = results
            print(f"  dps {dps:5d}  ->  speedup {m0/max(m1,1):.3f}x, "
                  f"anchors agree to {agree_digits(v0, v1)} digits")


if __name__ == "__main__":
    main()
