"""exp-057 candidate: remove the cancellation that IS the calibration loss.

isuperf walks y toward the fixed point L until |y - L| <= isuperfr, then
evaluates subst(fsl, x, y - L). At dps 400 isuperfr = 9.696e-20, so y and L are
both O(1) while their difference is at scale 1e-19: an absolute 1e-precis error
in y re-enters as a RELATIVE 1e-(precis-19) error in the argument. Everything
downstream inherits it.

Measured (this repo, base e, dps 400, precis = 404):

    seriesprecision   isuperfr    digits lost vs a ps=48 build
             21       9.70e-20            21.4
             32       3.80e-13            12.1

i.e. the loss is log10(1/isuperfr) = precis/21, and the engine delivers 383.8
true digits of precis = 404 end to end. The engine's accuracy IS isuperf's
accuracy, and the published calibration law ("dps - 24") is this cancellation
seen from outside.

The fix is a change of coordinate, not of algorithm: iterate u = y - L directly.
For x2mode == 0, with fs(z) = exp(z)-1+k, finv(z) = log(z-k+1) and
lambda1 = L-k+1 (line 235), the fixed point gives exp(L) = lambda1, hence

    fs(L+u)   - L = lambda1 * expm1(u)
    finv(L+u) - L = log1p(u / lambda1)

both exact, and the subtraction never happens. The same holds at the lower
fixed point with lambda2 = L2-k+1 (line 236). PARI 2.17.3 has complex-capable
log1p and expm1 at full relative accuracy (verified).

The initial u = z - L is harmless: z sits O(1) away from L (the walk is what
brings it close), so there is no cancellation there.

x2mode (the sub-eta / Shell-Thron regime) has a different finv with a sqrt
branch, so it keeps the original path untouched.

WARNING for the A/B: at FIXED dps this can look like a regression, because the
contour residual no longer stalls where it used to and the loop runs further.
Score it at fixed TRUE DIGITS (research/tools/digits_vs_reference.py), not at
fixed dps.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = Path(__file__).resolve().parent / "_exp057.gp"

PAT_ISUPERF = """isuperf(z) = {
  local(n,y,y1);
  y=z;
  n=0;
  if (repelling,
    while (abs(y-L)>isuperfr,
      y1=finv(y,L);
      if (abs(y1+2*Pi*I-L)<abs(y1-L), y=y1+2*Pi*I, y=y1);
      n++;
    );
  ,
    while (abs(y-L)>isuperfr,
      y=fs(y);
      n--;
    );
  );
  y = subst(fsl,x,y-L);"""

REP_ISUPERF = """isuperf(z) = {
  local(n,y,y1,u,u1);
  /* exp-057: iterate the OFFSET u = y - L instead of y. The old form walked y
     to within isuperfr = 10^(-precis/21) of L and only then formed y-L, which
     cancels log10(1/isuperfr) digits (measured: 21.4 at dps 400) -- that
     cancellation is the engine's whole calibration loss. In u-coordinates it
     never happens, using exp(L) = lambda1 = L-k+1:
         fs(L+u)   - L = lambda1*expm1(u)
         finv(L+u) - L = log1p(u/lambda1)
     x2mode has a different finv (sqrt branch) and keeps the original path. */
  if (x2mode==0,
    u = z - L;
    n = 0;
    if (repelling,
      while (abs(u)>isuperfr,
        u1 = log1p(u/lambda1);
        if (abs(u1+2*Pi*I)<abs(u1), u = u1+2*Pi*I, u = u1);
        n++;
      );
    ,
      while (abs(u)>isuperfr,
        u = lambda1*expm1(u);
        n--;
      );
    );
    y = subst(fsl,x,u);
  ,
  y=z;
  n=0;
  if (repelling,
    while (abs(y-L)>isuperfr,
      y1=finv(y,L);
      if (abs(y1+2*Pi*I-L)<abs(y1-L), y=y1+2*Pi*I, y=y1);
      n++;
    );
  ,
    while (abs(y-L)>isuperfr,
      y=fs(y);
      n--;
    );
  );
  y = subst(fsl,x,y-L);
  );"""

PAT_ISUPERF2 = """isuperf2(z) = {
  local(n,y);
  y=z;
  n=0;
  while (abs(y-L2)>isuperfr2,
    y=finv(y,L2);
    n++;
  );
  y = subst(fsl2,x,y-L2);"""

REP_ISUPERF2 = """isuperf2(z) = {
  local(n,y,u);
  /* exp-057: same offset iteration at the lower fixed point, lambda2 = L2-k+1 */
  if (x2mode==0,
    u = z - L2;
    n = 0;
    while (abs(u)>isuperfr2,
      u = log1p(u/lambda2);
      n++;
    );
    y = subst(fsl2,x,u);
  ,
  y=z;
  n=0;
  while (abs(y-L2)>isuperfr2,
    y=finv(y,L2);
    n++;
  );
  y = subst(fsl2,x,y-L2);
  );"""

PATCHES = [(PAT_ISUPERF, REP_ISUPERF), (PAT_ISUPERF2, REP_ISUPERF2)]


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
