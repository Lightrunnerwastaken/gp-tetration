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
        var: str = "ctrmul", timeout: float | None = None) -> dict:
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
        # REFERENCE-FREE guard. The anchor column below only works when a
        # ctrmul=1 run is in the same scan; this catches a bad setting on its
        # own. It has to be a ROUNDTRIP: the Abel equation A(f(z)) = A(z)+1 is
        # invariant under A -> A + theta(A) for any 1-periodic theta, so a run
        # that converged to the wrong theta still satisfies it (measured, base
        # 10 at ctrmul 0.85: Abel residual 0.E-96 while the value carried ZERO
        # correct digits). The engine's own `re` is a consistency signal, not a
        # correctness signal, for exactly the same reason.
        'print("SCAN-RT ", abs(sexp(slog(2.5)) - 2.5));',
        "quit;",
    ]) + "\n"
    t0 = time.time()
    # A degenerate setting does not always crash: with the exp-063 thsamples
    # clamp in place, ctrmul <= 0.5 at dps 300 keeps `re` creeping upward by
    # microscopic amounts while nlim extends itself (fork line ~1784), so the
    # run never terminates. Measured: 1855 s CPU against a 72 s baseline at
    # only 17 MB resident -- no large ct is being built, the loop just spins.
    # Without a cap one bad point hangs the whole scan and, worse, silently
    # discards every row that already completed.
    try:
        p = subprocess.run([GP, "-q", "-f"], input=script, capture_output=True,
                           text=True, encoding="utf-8", errors="replace",
                           cwd=str(REPO), timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"mul": mul, "timeout": round(time.time() - t0, 1)}
    ms = _re.search(r"SCAN-MS\s+(\d+)", p.stdout)
    tm = _re.search(r"SCAN-TERMS\s+(\d+)", p.stdout)
    vl = _re.search(r"SCAN-VAL\s+([-\d.eE]+)", p.stdout)
    it = _re.search(r"SCAN-ITER \[[^,]*, *(\d+)", p.stdout)
    # PARI prints "3.74 E-96" with a space before the exponent -- match it.
    rt = _re.search(r"SCAN-RT\s+([-\d.]+(?:\s*E-?\d+)?)", p.stdout)
    if not (ms and vl):
        return {"mul": mul, "error": p.stdout[-700:] + p.stderr[-400:]}
    return {"mul": mul, "ms": int(ms.group(1)), "deg": int(tm.group(1)) if tm else -1,
            "iter": int(it.group(1)) if it else -1,
            "val": vl.group(1), "wall": round(time.time() - t0, 1),
            "rt": rt.group(1).replace(" ", "") if rt else None}


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
    ap.add_argument("--timeout", type=float, default=900.0,
                    help="per-run wall-clock cap in seconds (0 = no cap)")
    args = ap.parse_args()
    eng = Path(args.engine)
    ref = None
    print(f"{args.var:>8s} {'init s':>9s} {'deg(ct)':>9s} {'iters':>7s} "
          f"{'digits vs 1':>13s} {'roundtrip':>12s}", flush=True)
    for mul in args.mul:
        r = run(mul, args.dps, args.base, eng, args.set, args.var,
                args.timeout or None)
        if "timeout" in r:
            print(f"{mul:>8s} {r['timeout']:>9.1f}  TIMEOUT (no convergence "
                  f"within cap; not a speed result)", flush=True)
            continue
        if "error" in r:
            print(f"{mul:>8s}  ERROR: {r['error'][:200]}", flush=True)
            continue
        if ref is None:
            ref = r["val"]
            d = "(reference)"
        else:
            d = str(agree(ref, r["val"]))
        rt, flag = r.get("rt"), ""
        if rt is not None:
            try:
                v = float(rt)
                flag = "  <-- BROKEN" if v > 1e-20 else ""
                rt = f"{v:.2e}"
            except ValueError:
                flag = "  <-- unparsed"
        print(f"{mul:>8s} {r['ms']/1000:9.2f} {r['deg']:9d} {r['iter']:7d} "
              f"{d:>13s} {rt or '-':>12s}{flag}", flush=True)


if __name__ == "__main__":
    main()
