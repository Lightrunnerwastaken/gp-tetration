"""Scan the sampling radius knob ctrmul -- no fork mutation needed.

loop() has `if (ctrmul != 1, ctr = ctr*ctrmul)`, so the sampling radius can be
scanned from outside. The cost law is T ~ I(p) * N(p)^2 * M(q) with

    N   = terms_per_digit * re,   terms_per_digit = 1 / log10(1/ctr)
    I   = re / rate(ctr)

so shrinking ctr shrinks N quadratically but slows the per-iteration rate.
exp-023 found ctr*0.9 optimal, but scanned only near 1 and at dps 150-400.
Since N enters squared and terms_per_digit falls fast (0.72 -> 7.0 terms per
digit, 0.5 -> 3.3), the product optimum may sit much lower -- worth checking
before any expensive evaluation-primitive work.

Reports init time AND the anchor, so a "win" that merely converged less far is
visible: the anchor is compared digit-by-digit against the ctrmul=1 run.

    python research/tools/ctr_scan.py --dps 200 --mul 1 0.9 0.8 0.7 0.6 0.5
"""
from __future__ import annotations

import argparse
import re as _re
import subprocess
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
FORK = REPO / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
GP = r"C:\Program Files (x86)\Pari64-2-17-3\gp.exe"


def run(mul: str, dps: int, base: str, engine: Path, extra: list[str] | None = None,
        var: str = "ctrmul") -> dict:
    script = "\n".join([
        "default(parisizemax, 8589934592);",
        f"default(realprecision, {dps});",
        f'read("{engine.as_posix()}");',
        f"default(realprecision, {dps});",
        "quietmode=1;",
        f"{var} = {mul};",
        *[f"{e};" for e in (extra or [])],
        f"sbase = {base};",
        "gettime();",
        "sres = sexpinit(sbase,0,0,0);",
        'print("SCAN-MS ", gettime());',
        'print("SCAN-ITER ", sres);',
        'print("SCAN-TERMS ", poldegree(ct));',
        'print("SCAN-VAL ", sexp(0.5));',
        "quit;",
    ]) + "\n"
    t0 = time.time()
    p = subprocess.run([GP, "-q", "-f"], input=script, capture_output=True,
                       text=True, encoding="utf-8", errors="replace", cwd=str(REPO))
    ms = _re.search(r"SCAN-MS\s+(\d+)", p.stdout)
    tm = _re.search(r"SCAN-TERMS\s+(\d+)", p.stdout)
    vl = _re.search(r"SCAN-VAL\s+([-\d.eE]+)", p.stdout)
    it = _re.search(r"SCAN-ITER \[[^,]*, *(\d+)", p.stdout)
    if not (ms and vl):
        return {"mul": mul, "error": p.stdout[-700:] + p.stderr[-400:]}
    return {"mul": mul, "ms": int(ms.group(1)), "deg": int(tm.group(1)) if tm else -1,
            "iter": int(it.group(1)) if it else -1,
            "val": vl.group(1), "wall": round(time.time() - t0, 1)}


def agree(a: str, b: str) -> int:
    n = 0
    for ca, cb in zip(a, b):
        if ca != cb:
            break
        if ca.isdigit():
            n += 1
    return n


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dps", type=int, default=200)
    ap.add_argument("--base", default="exp(1)")
    ap.add_argument("--mul", nargs="+", default=["1", "0.9", "0.8", "0.7", "0.6"])
    ap.add_argument("--engine", default=str(FORK))
    ap.add_argument("--var", default="ctrmul",
                    help="global to scan (ctrmul, irmul, ...)")
    ap.add_argument("--set", action="append", default=[],
                    help="extra GP assignment, e.g. --set irmul=0.8")
    args = ap.parse_args()
    eng = Path(args.engine)
    ref = None
    print(f"{args.var:>8s} {'init s':>9s} {'deg(ct)':>9s} {'iters':>7s} {'digits vs 1':>13s}")
    for mul in args.mul:
        r = run(mul, args.dps, args.base, eng, args.set, args.var)
        if "error" in r:
            print(f"{mul:>8s}  ERROR: {r['error'][:200]}")
            continue
        if ref is None:
            ref = r["val"]
            d = "(reference)"
        else:
            d = str(agree(ref, r["val"]))
        print(f"{mul:>8s} {r['ms']/1000:9.2f} {r['deg']:9d} {r['iter']:7d} {d:>13s}")


if __name__ == "__main__":
    main()
