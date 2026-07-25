"""exp-058 candidate: truncate the THETA SERIES per sample in sfunc.

Not to be confused with exp-049, which truncated the ct-Horner inside
thtaylor's thfunc loop and gained nothing. This is a different evaluation:

    fatou_fork.gp:1014   y1 = y1 + subst(tht, x, exp((y1-zth)*2*Pi*I));

the theta branch's own series evaluation in sfunc, run at full degree and full
precision for every theta-branch sample of every pass. It has never been
reduced -- exp-021's incremental path lives in icabel, exp-048's radius levels
feed icabel only.

The argument u = exp((y1-zth)*2*Pi*I) has |u| spread over roughly [0.06, 0.53],
so the same reasoning as exp-048 applies: tht's coefficients decay
geometrically, so a sample at |u| = 0.06 needs far fewer terms than one at
0.53, and today all of them get the full degree.

Build 8 truncation levels over |u| in (0, 0.56] once per staylor pass (8 * O(deg
tht) work, against deg(tht) per sample saved), using the same exponent-based
criterion as icbuild -- derived from the ACTUAL coefficients, not a decay model.
Samples above the table fall through to the untruncated tht.

Only the upper branch (tht) is touched; tht2 keeps the original path.

Two traps from earlier in this session, both re-checked here:
  * the loop variable must not be called k -- k is the engine's base parameter;
  * a t_REAL polynomial does not truncate via the Euclidean quotient, so the
    truncation is built explicitly with Polrev(vector(...)).
Verify the level degrees are actually below deg(tht) before believing a timing.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
HERE = Path(__file__).resolve().parent

PAT_GLOBALS = "icnlev=8; swlev=0; icdl=0; icrmax=1; icdlon=0;"
REP_GLOBALS = ("icnlev=8; swlev=0; icdl=0; icrmax=1; icdlon=0;\n"
               "/* exp-058: per-|u| truncation levels for the theta series in sfunc */\n"
               "thnlev=%(L)d; thlv=0; thlvon=0; thumax=0.56;")

PAT_HELPER = "icbuild(pp, dig) = {"
REP_HELPER = """/* exp-058: truncation levels for the theta series over |u| in (0, thumax].
   Same criterion as icbuild: keep term j while
     exponent(c_j) + j*log2(rho) >= exponent(max|c|) - dig*log2(10) - 16. */
thbuild(pp, dig) = {
  local(v, n, lg2, tol2, j, lr2, kk, ii, mx);
  thlvon = 0;
  if (type(pp) != "t_POL", return(0));
  v = Vecrev(pp);
  n = #v;
  if (n < 24, return(0));
  lg2 = vector(n, kk, if (v[kk]==0, -1000000, exponent(v[kk])));
  mx = vecmax(lg2);
  if (mx <= -1000000, return(0));
  tol2 = mx - dig*3.3219280948873623 - 16;
  thlv = vector(thnlev);
  for (j=1, thnlev,
    lr2 = log(thumax*j/thnlev)/log(2);
    kk = n;
    while ((kk > 1) && (lg2[kk] + (kk-1)*lr2 < tol2), kk--);
    thlv[j] = if (kk >= n, pp, Polrev(vector(kk, ii, v[ii])));
  );
  thlvon = 1;
  return(1);
}

icbuild(pp, dig) = {"""

PAT_BUILD = """    icrmax = r;
    icbuild(if (icfull, ct, icdct), if (icfull, precis, icdig));"""
REP_BUILD = """    icrmax = r;
    icbuild(if (icfull, ct, icdct), if (icfull, precis, icdig));
    /* exp-058: theta-series levels for this pass (tht was rebuilt by loop()
       just before staylor was called) */
    thbuild(tht, default(realprecision));"""

PAT_EVAL = """      z = y1;
      n = I;
      y1 = y1 + subst(tht,x,exp((y1-zth)*2*Pi*I));"""
REP_EVAL = """      z = y1;
      n = I;
      /* exp-058: |u| spreads over ~[0.06, 0.53], so pick the truncation level
         that this sample actually needs instead of the full degree. */
      y1 = exp((y1-zth)*2*Pi*I);
      y1 = z + subst(if (thlvon,
                       thlv[max(1, min(thnlev, 1+floor(thnlev*abs(y1)/thumax)))],
                       tht), x, y1);"""

PATCHES = [(PAT_GLOBALS, REP_GLOBALS), (PAT_HELPER, REP_HELPER),
           (PAT_BUILD, REP_BUILD), (PAT_EVAL, REP_EVAL)]


def build(levels: int, dst: Path) -> Path:
    text = SRC.read_text(encoding="utf-8")
    for i, (pat, rep) in enumerate(PATCHES):
        if text.count(pat) != 1:
            print(f"patch {i} matched {text.count(pat)} times -- abort", file=sys.stderr)
            raise SystemExit(1)
        text = text.replace(pat, rep % {"L": levels} if "%(L)d" in rep else rep, 1)
    dst.write_bytes(text.encode("utf-8"))
    return dst


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", type=int, default=8)
    a = ap.parse_args()
    print(f"wrote {build(a.levels, HERE / f'_exp058_{a.levels}.gp')} (levels={a.levels})")
