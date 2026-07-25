"""exp-041 candidate: evaluate the theta branch's y2 at reduced precision.

In sfunc() the theta branch computes y2 = icabel(zc) -- a full O(ctsamples)
Horner at icdig digits -- and uses it for THREE DISCRETE DECISIONS ONLY:

    if (imag(y2)>0, ...)                          branch (upper/lower theta)
    while (abs(y1-y2) > abs(y1-y2 +- Period), ...) period representative

Measured margins (research/tools/branchprobe.py, base e, dps 100/200/300):

    |imag(y2)|                >= 0.795     (branch decision)
    relative period-tie margin >= 0.9998    (period decision)
    period shifts performed    == 0         at every dps

So the decisions are settled by O(1) margins and survive any precision above
a few digits. This mutation gives icabel an optional precision argument and
lets that one call site pass it; every caller that needs the VALUE (the
direct branch, thest, abel) is untouched.

Theta-branch share: 33.6% of all samples at dps 100/200/300.

Writes research/tools/_exp041.gp -- no fork mutation until it is measured.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = Path(__file__).resolve().parent / "_exp041.gp"

PAT_GLOBAL = "icct=0; icvals=0; icA=0; ickey=0; icdct=0; icdig=60; icfull=1;"
REP_GLOBAL = (
    "icct=0; icvals=0; icA=0; ickey=0; icdct=0; icdig=60; icfull=1;\n"
    "/* exp-041: digits used for the theta branch's decision-only icabel call */\n"
    "y2dig=%d;"
)

PAT_ICABEL = """icabel(zz) = {
  local(h, A);
  if (swon && swidx>0,
    if (icfull,
      A = rlnlm*(log(I*(zz-L))-Pi*I/2) + rlnlm2*(log(-I*(zz-L2))+Pi*I/2) + sfunczero;
      h = subst(ct, x, (zz-circc));
      icA[swidx] = A;
      icvals[swidx] = h;
    ,
      h = icvals[swidx] + subst(icdct, x, precision(zz-circc, icdig));
      icvals[swidx] = h;
      A = icA[swidx];
    );"""

REP_ICABEL = """icabel(zz, pdig=0) = {
  local(h, A, zq);
  if (swon && swidx>0,
    if (icfull,
      A = rlnlm*(log(I*(zz-L))-Pi*I/2) + rlnlm2*(log(-I*(zz-L2))+Pi*I/2) + sfunczero;
      zq = if (pdig, precision(zz-circc, pdig), zz-circc);
      h = subst(ct, x, zq);
      icA[swidx] = A;
      icvals[swidx] = h;
    ,
      zq = precision(zz-circc, if (pdig, min(pdig, icdig), icdig));
      h = icvals[swidx] + subst(icdct, x, zq);
      icvals[swidx] = h;
      A = icA[swidx];
    );"""

PAT_CALL = "    y2 = icabel(zc); /* exp-015: lazy — only the theta path needs it */"
REP_CALL = ("    /* exp-015 lazy; exp-041: decisions only -> reduced precision */\n"
            "    y2 = icabel(zc, y2dig);")


def build(y2dig: int, dst: Path = DST) -> Path:
    text = SRC.read_text(encoding="utf-8")
    for pat, rep in ((PAT_GLOBAL, REP_GLOBAL % y2dig),
                     (PAT_ICABEL, REP_ICABEL),
                     (PAT_CALL, REP_CALL)):
        if text.count(pat) != 1:
            print(f"pattern matched {text.count(pat)} times -- abort", file=sys.stderr)
            raise SystemExit(1)
        text = text.replace(pat, rep, 1)
    dst.write_bytes(text.encode("utf-8"))
    return dst


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--y2dig", type=int, default=30)
    ap.add_argument("--out", default=str(DST))
    a = ap.parse_args()
    p = build(a.y2dig, Path(a.out))
    print(f"wrote {p} (y2dig={a.y2dig})")
