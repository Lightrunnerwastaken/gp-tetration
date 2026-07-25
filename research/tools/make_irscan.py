"""Scan knob for the inner (direct-vs-theta) radius: adds a global `irmul`.

Geometry as actually implemented:

  initsch()  sets ir = 7/10, ctr = 4/5, and  ircircr = ir*ctr*circr = 0.56*circr
  loop()     then does  if (efam, ctr = ctr*81/100)  ->  ctr = 0.648
             ircircr is NOT recomputed, so the sampling radius r = ctr*circr
             and the branch radius ircircr are already decoupled. That is not
             an oversight: z0h/z0l are placed on the ircircr circle inside
             initsch, so recomputing ircircr afterwards would decouple it from
             the theta anchor instead. Scanning it down was measured and is
             worse -- irmul 0.85 collapses the run to 92 digits.

  sfunc:     |z - circc| <  ircircr  ->  DIRECT branch  (evaluates ct: O(N))
             |z - circc| >= ircircr  ->  THETA branch   (since exp-041b: no ct
                                          at all, only the degree-thsamples
                                          theta series, ~14x smaller)

Two consequences that the ctrmul scan exposed:

  * the sampling radius cannot fall below ircircr/circr = 0.56. **This floor
    moved with exp-045 and the margin is now thin.** With the exp-023 factor
    9/10 the radius sat at ctr = 0.72, i.e. 29% clear of the floor (ctrmul
    could go to 0.778). exp-045 doubled the factor to 81/100, so ctr = 0.648 --
    only 16% clear, and the ctrmul floor is now **0.864**, not 0.778. A scan
    that reuses the old headroom will walk off the cliff: below it the sampling
    circle lies INSIDE the branch radius, the theta fit degenerates and
    thsamples goes negative. Measured at the old radius: 0.8 runs, 0.7 dies.
  * lowering ircircr helps twice: more samples move to the (cheap) theta
    branch, AND the floor under ctrmul drops, allowing a smaller grid.

The risk is accuracy: ircircr is where the theta/Schroeder representation is
deemed valid. ir is auto-selected as 7/10 or 15/16 and matrix_ir even uses
24/25 -- i.e. the author only ever went LARGER. Going smaller is uncharted, so
the anchor agreement against irmul=1 is the acceptance test, not the clock.

Writes research/tools/_irscan.gp (scan instrument, not a keep candidate).
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = Path(__file__).resolve().parent / "_irscan.gp"

PATCHES = [
    ("ctrmul = 1;  /* exp-023 scan knob: scales the theta sampling radius */",
     "ctrmul = 1;  /* exp-023 scan knob: scales the theta sampling radius */\n"
     "irmul = 1;   /* scan knob: scales the direct-vs-theta branch radius */"),
    ("  if (ctrmul != 1, ctr = ctr*ctrmul);  /* scan knob, default 1 */",
     "  if (ctrmul != 1, ctr = ctr*ctrmul);  /* scan knob, default 1 */\n"
     "  if (irmul != 1, ircircr = ircircr*irmul);  /* scan knob, default 1 */"),
]


def main() -> int:
    text = SRC.read_text(encoding="utf-8")
    for i, (pat, rep) in enumerate(PATCHES):
        if text.count(pat) != 1:
            print(f"patch {i} matched {text.count(pat)} times -- abort", file=sys.stderr)
            return 1
        text = text.replace(pat, rep, 1)
    DST.write_bytes(text.encode("utf-8"))
    print(f"wrote {DST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
