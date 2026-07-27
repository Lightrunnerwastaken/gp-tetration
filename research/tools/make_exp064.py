"""exp-064: icbuild scans the coefficient array once per radius level; one sweep suffices.

exp-048 keeps `icnlev` truncated copies of the difference polynomial, one per
radius bucket, so that a sample at radius rho only evaluates the degree it
actually needs. Building them costs a scan per level:

    for (j=1, icnlev,
      lr2 = log(icrmax*j/icnlev)/log(2);
      kk = n;                                    <-- reset every level
      while ((kk > 1) && (lg2[kk] + (kk-1)*lr2 < tol2), kk--);
      ...)

so the array is walked from the top icnlev = 8 times per staylor pass.

**The cutoffs are provably monotone in the level.** lr2 is increasing in j
(icrmax*j/icnlev increases), hence lg2[kk] + (kk-1)*lr2 is increasing in j for
every kk > 1, hence the while-condition is harder to satisfy at larger j, hence
the retained length kk_j is non-decreasing in j. Sweeping j from icnlev DOWN to
1 and carrying kk therefore lands on exactly the same cutoffs: at level j the
carried kk equals kk_{j+1} >= kk_j, and the loop only decrements, so it walks
precisely the indices the fresh scan would have walked. One pass, same result.

Measured share of the run: 3.65% (doubling probe on the icbuild call site).

The second-order effect is the larger one. exp-048 rejected 16 and 32 radius
levels because of the build cost -- and this scan IS the build cost. With an
O(N) build rather than O(icnlev*N), icnlev becomes worth re-scanning, and
exp-048's own table shows the outer band (382 of 1280 points, still at full
degree) carries most of the remaining Horner work.

Acceptance test is bit-identity, not the clock: this changes no arithmetic.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = Path(__file__).resolve().parent / "_exp064.gp"

PAT = """  icdl = vector(icnlev);
  for (j=1, icnlev,
    lr2 = log(icrmax*j/icnlev)/log(2);
    kk = n;
    while ((kk > 1) && (lg2[kk] + (kk-1)*lr2 < tol2), kk--);"""

REP = """  icdl = vector(icnlev);
  /* exp-064: the cutoffs are monotone in the level. lr2 increases with j, so
     lg2[kk] + (kk-1)*lr2 increases with j, so the while-condition is harder to
     satisfy and the retained length kk_j is non-decreasing in j. Sweeping j
     DOWNWARD and carrying kk therefore reproduces the cutoffs exactly, in ONE
     pass over the array instead of icnlev passes. Bit-identical by
     construction; measured 3.65% of the run. */
  kk = n;
  forstep (j=icnlev, 1, -1,
    lr2 = log(icrmax*j/icnlev)/log(2);
    while ((kk > 1) && (lg2[kk] + (kk-1)*lr2 < tol2), kk--);"""


def main() -> int:
    text = SRC.read_text(encoding="utf-8")
    if text.count(PAT) != 1:
        print(f"anchor matched {text.count(PAT)} times -- abort", file=sys.stderr)
        return 1
    DST.write_bytes(text.replace(PAT, REP, 1).encode("utf-8"))
    print(f"wrote {DST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
