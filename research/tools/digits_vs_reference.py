"""Measure how many TRUE digits an engine delivers, against the proven references.

Comparing two engines against each other only shows where they diverge, not
which one is right. research/reference/values.json holds error-vector-proven
values (sexp_e(0.5) to 972 digits), so the honest question -- "did this
mutation cost real digits?" -- is answered here, not by the gate (whose e-group
thresholds are a low legacy of the v2 correction) and not by mutual agreement.

    python research/tools/digits_vs_reference.py --dps 300 \
        --engine src/fatou_backend/vendor/fatou_fork.gp \
        --engine research/tools/_exp044_1_60.gp
"""
from __future__ import annotations

import argparse
import json
import re as _re
import subprocess
import time
from pathlib import Path

import mpmath as mp

REPO = Path(__file__).resolve().parents[2]
GP = r"C:\Program Files (x86)\Pari64-2-17-3\gp.exe"
REF = REPO / "research" / "reference" / "values.json"


def reference(key: str) -> mp.mpf:
    payload = json.loads(REF.read_text(encoding="utf-8"))
    e = payload["values"][key]
    return mp.mpf(e["real"])


def run(engine: Path, dps: int, base: str, arg: str, extra: list[str]) -> tuple[str, int]:
    script = "\n".join([
        "default(parisizemax, 8589934592);",
        f"default(realprecision, {dps});",
        f'read("{engine.as_posix()}");',
        f"default(realprecision, {dps});",
        "quietmode=1;",
        *[f"{e};" for e in extra],
        f"rbase = {base};",
        "gettime();",
        "sexpinit(rbase,0,0,0);",
        'print("RMS ", gettime());',
        f'print("RVAL ", sexp({arg}));',
        "quit;",
    ]) + "\n"
    p = subprocess.run([GP, "-q", "-f"], input=script, capture_output=True,
                       text=True, encoding="utf-8", errors="replace", cwd=str(REPO))
    v = _re.search(r"RVAL\s+([-\d.eE ]+)", p.stdout)
    m = _re.search(r"RMS\s+(\d+)", p.stdout)
    if not v:
        raise SystemExit(p.stdout[-1200:] + p.stderr[-500:])
    return v.group(1).strip().replace(" E", "e").replace(" ", ""), int(m.group(1)) if m else -1


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", action="append", required=True)
    ap.add_argument("--dps", type=int, nargs="+", default=[300])
    ap.add_argument("--base", default="exp(1)")
    ap.add_argument("--arg", default="0.5")
    ap.add_argument("--key", default="sexp|e|0.5|1000")
    ap.add_argument("--set", action="append", default=[])
    args = ap.parse_args()

    mp.mp.dps = 1100
    ref = reference(args.key)
    print(f"reference {args.key} loaded ({mp.mp.dps} dps working)")
    for dps in args.dps:
        print(f"\n-- dps {dps} --")
        for e in args.engine:
            eng = Path(e) if Path(e).is_absolute() else REPO / e
            t0 = time.time()
            val, ms = run(eng, dps, args.base, args.arg, args.set)
            got = mp.mpf(val)
            err = abs(got - ref)
            digits = float(-mp.log10(err / max(abs(ref), mp.mpf(1)))) if err else float(dps)
            print(f"  {eng.name:26s} {ms/1000:8.2f}s   TRUE digits {digits:8.1f}"
                  f"   (wall {time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
