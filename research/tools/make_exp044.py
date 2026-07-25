"""exp-044 candidate: relax exp-008's precision ladder for the e-family.

loop() runs early iterations at reduced working precision (exp-008):

    default(realprecision, max(48, min(precis, 2*floor(re) + 60)));

The 2x factor is there for the FACTOR-2 bases: exp-006 used re+40 and "cut
real digits of the fast bases" (b=2, b=10 carry up to ~2*re true digits in the
state). The e-family does not have that property -- journal 2026-07-10 records
"sexp error == contour error, factor 1" -- so for efam the doubling buys
nothing and costs arithmetic precision in exactly the cost-dominant middle of
the run (r ~ 0.69p, where 2*re+60 has already saturated at precis while
re+60 would still be ~0.31p below it).

Expected: no effect at all below re ~ (precis-60)/2, then a growing win with
depth. dps 1020: at re=690 the working precision would be 750 instead of 1020.

DANGER: this is exactly the class of change that silently costs digits, and
the gate's e-group thresholds are LOW (58.9 / 59.2, a legacy of the v2
reference correction). Do NOT trust the gate alone here -- check
min_correct_digits in the benchmark (cap 190 for the dps-200 cases) and the
anchor agreement against the unmutated engine at dps 300.

    python research/tools/make_exp044.py --mult 1
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
HERE = Path(__file__).resolve().parent

PAT = "    default(realprecision, max(48, min(precis, 2*floor(re) + 60)));"
REP = """    /* exp-044: the 2x guard protects the FACTOR-2 bases, whose state
       carries up to 2*re true digits (exp-006's re+40 cut them). The
       e-family runs at factor 1.

       The point is not the constant but the SHAPE. icdig (and exdig, thdig)
       are computed as realprecision - re + 40, so with a factor m the
       increment precision is (m-1)*re + margin + 40: linear in re for m=2,
       CONSTANT for m=1. The cost sum sum_re N(re)^2 * M(icdig) therefore
       drops from ~p^3.9 to ~p^3.05 in the model.

       Accuracy constraint (this is what bounds the margin): within one grid
       stretch the incremental accumulator icvals only carries
       m*re_start + margin + 40 digits, and it must still cover the residual
       at the END of the stretch. So margin must exceed the largest stretch
       span in re. Measured spans at dps 200: 10..38 digits. */
    default(realprecision, max(48, min(precis, floor(%(m)s*re) + %(g)s)));"""


def build(mult: str, dst: Path, margin: str = "60") -> Path:
    text = SRC.read_text(encoding="utf-8")
    if text.count(PAT) != 1:
        print(f"pattern matched {text.count(PAT)} times -- abort", file=sys.stderr)
        raise SystemExit(1)
    expr = f"if(efam,{mult},2)"
    gexpr = f"if(efam,{margin},60)"
    dst.write_bytes(text.replace(PAT, REP % {"m": expr, "g": gexpr}, 1).encode("utf-8"))
    return dst


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mult", default="1", help="multiplier used for the e-family")
    ap.add_argument("--margin", default="60", help="additive guard for the e-family")
    a = ap.parse_args()
    tag = (a.mult + "_" + a.margin).replace(".", "").replace("/", "")
    p = build(a.mult, HERE / f"_exp044_{tag}.gp", a.margin)
    print(f"wrote {p} (efam ladder {a.mult}*re + {a.margin})")
