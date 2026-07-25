"""Run the section profiler (_prof.gp) and summarize where sexpinit's time goes.

    python research/tools/make_profile.py            # build _prof.gp first
    python research/tools/section_profile.py --dps 200 400

Reports per run: total ms, and the split into
  tsmp  staylor sampling loop (N points x degree-N Horner)
  tth   thtaylor (theta rebuild: its own Horner block + a quadratic DFT)
  text  staylor extraction (FFT / Bluestein)
plus the cold (grid-stretch start) vs warm share.

Windows timer granularity is ~15.6 ms, so single-iteration numbers are
quantized; only the sums over all iterations are meaningful.
"""
from __future__ import annotations

import argparse
import collections
import json
import re
import subprocess
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
ENGINE = REPO / "research" / "tools" / "_prof.gp"
GP = r"C:\Program Files (x86)\Pari64-2-17-3\gp.exe"


def run(dps: int, base: str, gp: str = GP) -> dict:
    script = "\n".join([
        "default(parisizemax, 8589934592);",
        f"default(realprecision, {dps});",
        f'read("{ENGINE.as_posix()}");',
        f"default(realprecision, {dps});",
        "quietmode=1;",
        f"pbase = {base};",
        "gettime();",
        "sexpinit(pbase,0,0,0);",
        'print("TOTAL-MS ", gettime());',
        'print("ANCHOR ", sexp(0.5));',
        "quit;",
    ]) + "\n"
    t0 = time.time()
    proc = subprocess.run([gp, "-q", "-f"], input=script, capture_output=True,
                          text=True, encoding="utf-8", errors="replace", cwd=str(REPO))
    rows = []
    for line in proc.stdout.splitlines():
        if line.startswith("PROF "):
            d = dict(re.findall(r"(\w+)=([-\d.eE]+)", line))
            rows.append({k: (float(v) if "." in v else int(v)) for k, v in d.items()})
    total = re.search(r"TOTAL-MS\s+(\d+)", proc.stdout)
    anchor = re.search(r"ANCHOR\s+([-\d.eE]+)", proc.stdout)
    if not rows:
        raise SystemExit(f"no PROF lines\n{proc.stdout[-1500:]}\n{proc.stderr[-800:]}")
    agg, cold, warm = collections.Counter(), collections.Counter(), collections.Counter()
    for x in rows:
        tgt = cold if x["icfull"] else warm
        for k in ("tth", "tsmp", "text", "tsta"):
            agg[k] += x[k]
            tgt[k] += x[k]
        tgt["n"] += 1
    tot_ms = int(total.group(1)) if total else 0
    return {
        "dps": dps, "base": base, "wall_s": round(time.time() - t0, 1),
        "total_ms": tot_ms, "iterations": len(rows),
        "final_ctsamples": rows[-1]["cts"], "final_thsamples": rows[-1]["ths"],
        "final_re": rows[-1]["re"],
        "sum": dict(agg), "cold": dict(cold), "warm": dict(warm),
        "pct": {k: round(100 * v / tot_ms, 1) for k, v in agg.items()} if tot_ms else {},
        "cold_pct": round(100 * cold["tsta"] / tot_ms, 1) if tot_ms else 0,
        "anchor": anchor.group(1)[:45] if anchor else None,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dps", type=int, nargs="+", default=[200])
    ap.add_argument("--base", default="exp(1)")
    ap.add_argument("--out")
    args = ap.parse_args()
    out = []
    for dps in args.dps:
        r = run(dps, args.base)
        out.append(r)
        print(json.dumps(r), flush=True)
    if args.out:
        Path(args.out).write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
