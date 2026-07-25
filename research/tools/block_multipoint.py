"""A5: is the product-tree guard loss per BLOCK? If so the exponent moves.

The earlier guard sweep in multipoint_bench.py was invalid: it built the arb
inputs once at a fixed precision and afterwards raised only the working
precision, so the frozen input ball radii pinned every result. Rebuilding the
inputs at each precision gives a clean law -- about 6 guard bits per point.

That opens the only exponent lever left. If the loss really is c*K bits for a
block of K points (rather than c*N for the whole set), then chunking the points
into N/K blocks caps the guard at g*K instead of g*N, and the per-iteration cost
becomes

    (N/K) * M_poly(N) * M(q + g*K)      instead of      N^2 * M(q)

With K chosen so that g*K ~ q -- i.e. K growing linearly in p -- the precision
never more than doubles and the operation count falls by a factor ~K/log N.
That is the Moroz mechanism without implementing Moroz, and it would take the
per-iteration cost from Theta(p^2) to Theta(p log p).

This script measures the two things that decide it:
  1. does a g*K guard fully recover a block of K points (accuracy), and
  2. what does the blocked evaluation actually cost against Horner (time).

Inputs are rebuilt at each precision -- that is the whole point.

    python research/tools/block_multipoint.py --state research/tools/state_300.txt \
        --qdigits 100 300 --blocks 128 256 512 1280 --guards 3 6
"""
from __future__ import annotations

import argparse
import math
import time
from pathlib import Path

from flint import acb_poly, ctx

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from multipoint_bench import State, build


def blocked_fast(poly: acb_poly, pts: list, K: int) -> list:
    out = []
    for i in range(0, len(pts), K):
        out.extend(poly.evaluate(pts[i:i + K], algorithm="fast"))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--state", required=True)
    ap.add_argument("--qdigits", type=int, nargs="+", default=[300])
    ap.add_argument("--blocks", type=int, nargs="+", default=[128, 256, 512, 1280])
    ap.add_argument("--guards", type=float, nargs="+", default=[3.0, 6.0])
    args = ap.parse_args()

    st = State(Path(args.state))
    print(f"state dps={st.dps}")

    for q in args.qdigits:
        qb = int(q * 3.3219)
        # truth: Horner with inputs built well above the target
        ctx.prec = qb + 200 * 3 + 256
        poly, pts, steps = build(st)
        N, M = len(poly.coeffs()), len(pts)
        t0 = time.time()
        truth = poly.evaluate(pts, algorithm="iter")
        t_truth = time.time() - t0

        # baseline: Horner at the working precision, inputs built there
        ctx.prec = qb
        poly, pts, _ = build(st)
        t0 = time.time()
        hv = poly.evaluate(pts, algorithm="iter")
        t_horner = time.time() - t0
        err_h = max(_rel(a, b) for a, b in zip(hv, truth))
        print(f"\n=== q = {q} digits ({qb} bits), deg = {N-1}, points = {M} ===")
        print(f"  truth (horner, +200 digits)   {t_truth:7.2f}s")
        print(f"  horner @ q                    {t_horner:7.2f}s   rel err {err_h:.2e}")
        print(f"  {'block K':>8s} {'guard':>7s} {'prec bits':>10s} {'time':>8s} "
              f"{'rel err':>11s} {'vs horner':>10s}")
        for g in args.guards:
            for K in args.blocks:
                if K > M:
                    continue
                guard = int(g * K)
                ctx.prec = qb + guard
                poly, pts, _ = build(st)      # <-- rebuilt AT this precision
                t0 = time.time()
                vals = blocked_fast(poly, pts, K)
                dt = time.time() - t0
                err = max(_rel(a, b) for a, b in zip(vals, truth))
                print(f"  {K:8d} {g:7.1f} {qb+guard:10d} {dt:7.2f}s "
                      f"{err:11.2e} {t_horner/dt:9.2f}x")


def _rel(a, b) -> float:
    try:
        d = float((a - b).abs_upper())
        m = float(b.abs_lower())
        return d / m if m > 0 else d
    except Exception:
        return float("inf")


if __name__ == "__main__":
    main()
