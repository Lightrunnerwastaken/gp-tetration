"""F1: benchmark evaluation primitives for the operation that dominates ~88%
of sexpinit -- evaluating the degree-N polynomial ct at ~N scattered points.

Runs on a REAL state dump (research/tools/dump_state.py), because the point
distribution is exactly what E1b tripped over: the walk endpoints spread over
the whole disk (|w| ~ 0.07 .. 0.96), so a plain product tree underflows by
~N*log2(1/median) bits. E1b used 1024 guard bits where >=1766 were needed at
N=2048, and it measured at dps 300, below any plausible crossover.

Methods
  horner   FLINT's iterative evaluation (same algorithm GP uses)  -- baseline
  fast     product-tree multipoint at working precision + guard bits
  shells   points grouped by |w| into geometric shells, each shell rescaled to
           |u| ~ 1 (polynomial rescaled the same way) before its own product
           tree -- removes the underflow structurally

Truth is `horner` at working precision + 200 digits. Reported per method:
wall time, max relative error vs truth, and the worst ball radius, so a fast
method that is merely *fast* cannot pass as a win.

    python research/tools/multipoint_bench.py --state research/tools/state_300.txt \
        --qdigits 100 200 --guard-mult 2.0
"""
from __future__ import annotations

import argparse
import math
import time
from pathlib import Path

from flint import acb, acb_poly, arb, ctx


# ---------------------------------------------------------------- state dump

def _num(tok: str) -> arb:
    """PARI decimal, possibly in its '1.23... E-41' form, as an exact arb."""
    return arb(tok.replace(" E", "e").replace(" ", ""))


class State:
    def __init__(self, path: Path):
        self.meta: dict[str, str] = {}
        self.ct: list[tuple[str, str]] = []
        self.sw: list[tuple[int, int, str, str]] = []
        self.th: list[tuple[str, str]] = []
        for line in path.read_text().splitlines():
            if not line or line.startswith("#"):
                continue
            k, _, rest = line.partition(" ")
            rest = rest.strip()
            if k == "ct":
                a, _, b = rest.partition("|")
                self.ct.append((a, b))
            elif k == "sw":
                n, valid, pair = rest.split(" ", 2)
                a, _, b = pair.partition("|")
                self.sw.append((int(n), int(valid), a, b))
            elif k == "th":
                a, _, b = rest.partition("|")
                self.th.append((a, b))
            else:
                self.meta[k] = rest
        self.dps = int(self.meta["dps"])
        self.circc = tuple(self.meta["circc"].split("|"))
        self.sampr = tuple(self.meta["sampr"].split("|"))


def _acb(re_s: str, im_s: str) -> acb:
    return acb(_num(re_s), _num(im_s))


def build(st: State) -> tuple[acb_poly, list[acb], list[int]]:
    poly = acb_poly([_acb(r, i) for r, i in st.ct])
    circc = _acb(*st.circc)
    pts, steps = [], []
    for n, valid, r, i in st.sw:
        if not valid:
            continue
        pts.append(_acb(r, i) - circc)          # ct is evaluated at (z - circc)
        steps.append(n)
    for r, i in st.th:
        z = _acb(r, i)
        if z.abs_upper() == 0:
            continue
        pts.append(z - circc)
        steps.append(9999)                       # theta batch marker
    return poly, pts, steps


# ---------------------------------------------------------------- primitives

def ev_horner(poly: acb_poly, pts: list[acb]) -> list[acb]:
    return poly.evaluate(pts, algorithm="iter")


def ev_fast(poly: acb_poly, pts: list[acb]) -> list[acb]:
    return poly.evaluate(pts, algorithm="fast")


SHELL_RATIO = 2.0


def ev_shells(poly: acb_poly, pts: list[acb], ratio: float = 0.0) -> list[acb]:
    ratio = ratio or SHELL_RATIO
    """Group by |w| into geometric shells, rescale each shell to |u| ~ 1."""
    mags = [float(p.abs_lower()) for p in pts]
    lo = min(m for m in mags if m > 0)
    idx_by_shell: dict[int, list[int]] = {}
    for i, m in enumerate(mags):
        s = 0 if m <= 0 else int(math.log(m / lo) / math.log(ratio))
        idx_by_shell.setdefault(s, []).append(i)
    out: list[acb] = [acb(0)] * len(pts)
    coeffs = poly.coeffs()
    for s, idx in idx_by_shell.items():
        rho = arb(lo * ratio ** (s + 0.5))
        # q(u) = p(rho*u): coefficients a_k * rho^k  (one O(N) pass per shell)
        pw = arb(1)
        qc = []
        for c in coeffs:
            qc.append(c * pw)
            pw = pw * rho
        q = acb_poly(qc)
        us = [pts[i] / rho for i in idx]
        vals = q.evaluate(us, algorithm="fast")
        for i, v in zip(idx, vals):
            out[i] = v
    return out


METHODS = {"horner": ev_horner, "fast": ev_fast, "shells": ev_shells}


# ---------------------------------------------------------------- comparison

def rel_width(v: acb) -> float:
    """Ball width relative to the value -- the honest 'how many digits survive'."""
    try:
        hi, lo = float(v.abs_upper()), float(v.abs_lower())
        return float("inf") if lo <= 0 else (hi - lo) / lo
    except Exception:
        return float("inf")


def rel_err(a: acb, b: acb) -> float:
    d = (a - b).abs_upper()
    m = b.abs_lower()
    try:
        return float(d) / float(m) if float(m) > 0 else float(d)
    except Exception:
        return float("inf")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--state", required=True)
    ap.add_argument("--qdigits", type=int, nargs="+", default=[200])
    ap.add_argument("--guard-mult", type=float, default=2.0,
                    help="guard bits = guard_mult * N (van der Hoeven's O(N))")
    ap.add_argument("--methods", nargs="+", default=["horner", "fast", "shells"])
    ap.add_argument("--limit-points", type=int, default=0)
    ap.add_argument("--shell-ratio", type=float, default=2.0)
    args = ap.parse_args()

    globals()["SHELL_RATIO"] = args.shell_ratio
    st = State(Path(args.state))
    ctx.prec = int((st.dps + 60) * 3.3219) + 64
    poly, pts, steps = build(st)
    if args.limit_points:
        pts, steps = pts[:args.limit_points], steps[:args.limit_points]
    N = len(poly.coeffs())
    mags = sorted(float(p.abs_lower()) for p in pts)
    med = mags[len(mags) // 2]
    arcs = sorted(set(steps))
    print(f"state dps={st.dps}  deg(ct)={N-1}  points={len(pts)}  "
          f"|w| in [{mags[0]:.4f}, {mags[-1]:.4f}] median {med:.4f}")
    print(f"distinct walk-step counts (arcs): {len(arcs)}  -> {arcs[:20]}")
    print(f"product-tree underflow at N={len(pts)}: "
          f"{len(pts)*math.log2(1/max(med,1e-9)):.0f} bits")

    for q in args.qdigits:
        base_bits = int(q * 3.3219)
        guard = int(args.guard_mult * len(pts))
        print(f"\n--- working precision {q} digits ({base_bits} bits), "
              f"guard {guard} bits ---")
        ctx.prec = base_bits + 200 * 3
        t0 = time.time()
        truth = ev_horner(poly, pts)
        t_truth = time.time() - t0
        print(f"  truth (horner @ +200 digits): {t_truth:8.2f}s")
        for name in args.methods:
            ctx.prec = base_bits + (guard if name in ("fast", "shells") else 0)
            fn = METHODS[name]
            t0 = time.time()
            vals = fn(poly, pts)
            dt = time.time() - t0
            err = max(rel_err(v, t) for v, t in zip(vals, truth))
            width = max(rel_width(v) for v in vals)
            print(f"  {name:8s} {dt:8.2f}s   max rel err {err:.3e}   "
                  f"max rel ball width {width:.3e}")


if __name__ == "__main__":
    main()
