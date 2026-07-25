"""exp-046 candidate: keep thsamples in range so smaller sampling radii can run.

The ctrmul scan (exp-045) dies at ctrmul <= 0.7 with

    *** in function thtaylor: ...terms=samples-1;t_est=vector(samples,i,0)
    *** vector: domain error in vector: dimension < 0

because loop() computes

    thsamples = floor(re*tht_re_mult + 6);
    ...
    thlog = abs(polcoeff(tht,m));
    if (thlog<>0, thlog = log(thlog)+0.69;
                  tht_re_mult = -m*1.05*log10/thlog);

If the theta coefficient at index m has not decayed below e^-0.69 ~ 0.50, then
thlog > 0 and tht_re_mult turns NEGATIVE, so thsamples goes negative and the
vector allocation fails. A negative multiplier means "the theta series needs
MORE terms", so the right response is to grow the grid, not to shrink it.

This is a robustness fix, not a speed change: at the shipped radius the branch
never fires (the scan runs clean at ctrmul 0.8..1.0), so the engine must be
bit-identical there. That is the acceptance criterion -- verify the anchor is
unchanged at ctrmul=1 before believing any speed number at a lower radius.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = Path(__file__).resolve().parent / "_exp046.gp"

PAT = """      thsamples=floor(re*tht_re_mult+6);"""
REP = """      /* exp-046: tht_re_mult < 0 means the theta coefficient at the probe
         index has NOT decayed, i.e. the grid is too small -- growing it is
         the correct response. The old form produced a negative thsamples and
         crashed thtaylor's vector allocation. Never fires at the shipped
         radius; only smaller sampling radii reach it. */
      if (tht_re_mult > 0,
        thsamples=floor(re*tht_re_mult+6)
      ,
        thsamples=max(16, 2*thsamples));"""


def main() -> int:
    text = SRC.read_text(encoding="utf-8")
    if text.count(PAT) != 1:
        print(f"pattern matched {text.count(PAT)} times -- abort", file=sys.stderr)
        return 1
    DST.write_bytes(text.replace(PAT, REP, 1).encode("utf-8"))
    print(f"wrote {DST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
