"""exp-048 candidate: evaluate ct only to the degree each point actually needs.

Both Horner blocks -- staylor's sampling loop (67.9% of the run at dps 400) and
thtaylor's thfunc block (21.4%) -- evaluate the FULL polynomial at every point.
But ct's coefficients decay like circr^-k, so at a point of radius rho the tail
beyond term n contributes ~(rho/circr)^n: a point at rho = 0.55 needs 2.7x
fewer terms than a point at the sampling radius 0.96, for the same absolute
accuracy.

Measured on the real dps-300 state dump (research/tools/dump_state.py), tail
below 1e-300:

    |w| bucket   points   terms needed        (full degree in use: 2560)
    [0.0,0.2)      100        272
    [0.2,0.4)      365        415
    [0.4,0.6)      234        660
    [0.6,0.8)      199       1000
    [0.8,1.0)      382       2561

    total Horner work 3,276,800 -> 1,377,921   =  2.38x less

Implementation: bucket the samples by |z - circc| into `icnlev` levels (the
level is a property of the cached walk endpoint, so it is computed once per
grid stretch alongside swz), and keep one truncated copy of the polynomial per
level, rebuilt once per pass. The truncation index is derived from the ACTUAL
coefficient magnitudes via exponent() -- cheap integer work, no logs -- so it
adapts to whatever the diff polynomial looks like instead of trusting a decay
model. Building the copies costs icnlev * O(N) cheap operations against the
O(N^2) they save.

Safety: the truncation threshold is the largest coefficient scaled down by the
requested digits, with a 16-bit margin, and the level radius used is the UPPER
edge of the bucket, so every point in a bucket is covered conservatively.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
HERE = Path(__file__).resolve().parent

PAT_GLOBALS = "icct=0; icvals=0; icA=0; ickey=0; icdct=0; icdig=60; icfull=1;"
REP_GLOBALS = """icct=0; icvals=0; icA=0; ickey=0; icdct=0; icdig=60; icfull=1;
/* exp-048: per-point degree truncation. ct's coefficients decay like
   circr^-k, so a sample at radius rho needs only ~dig/log10(circr/rho) terms.
   Samples are bucketed by radius (cached per grid stretch, like swz) and each
   level keeps a truncated copy of the polynomial, rebuilt once per pass. */
icnlev=%(L)d; swlev=0; icdl=0; icrmax=1; icdlon=0;"""

# build the per-level truncations; placed right before staylor()
PAT_FUNC = "staylor( w,r,samples) = {"
REP_FUNC = """/* exp-048: build one truncated copy of pp per radius level.
   Truncation index from the actual coefficient magnitudes via exponent()
   (binary exponents, no logs): keep term k while
       exponent(c_k) + k*log2(rho)  >=  exponent(max|c|) - dig*log2(10) - 16.
   The level radius is the UPPER edge of the bucket, so every sample assigned
   to that level is covered. */
icbuild(pp, dig) = {
  /* NOTE: the loop variable must NOT be called k -- k is the base parameter
     of this engine (fs(z) = e^z - 1 + k), and using it here silently defeats
     the truncation. */
  local(v, n, lg2, tol2, j, lr2, kk, ii, mx);
  icdlon = 0;
  if (type(pp) != "t_POL", return(0));
  v = Vecrev(pp);
  n = #v;
  if (n < 64, return(0));
  lg2 = vector(n, kk, if (v[kk]==0, -1000000, exponent(v[kk])));
  mx = vecmax(lg2);
  if (mx <= -1000000, return(0));
  tol2 = mx - dig*3.3219280948873623 - 16;
  icdl = vector(icnlev);
  for (j=1, icnlev,
    lr2 = log(icrmax*j/icnlev)/log(2);
    kk = n;
    while ((kk > 1) && (lg2[kk] + (kk-1)*lr2 < tol2), kk--);
    /* build the truncation explicitly from the coefficient vector: the
       polynomial Euclidean quotient does not truncate a t_REAL polynomial
       the way one expects. */
    icdl[j] = if (kk >= n, pp, Polrev(vector(kk, ii, v[ii])));
  );
  icdlon = 1;
  return(1);
}

staylor( w,r,samples) = {"""

PAT_ALLOC = """      swisf = vector(samples);
      swisfv = vector(samples);"""
REP_ALLOC = """      swisf = vector(samples);
      swisfv = vector(samples);
      swlev = vector(samples);   /* exp-048: per-sample radius level */"""

PAT_BUILD = """    ickey = [samples, w, r];
    icct = ct;
  );"""
REP_BUILD = """    ickey = [samples, w, r];
    icct = ct;
    /* exp-048: (re)build the per-level truncations for this pass */
    icrmax = r;
    icbuild(if (icfull, ct, icdct), if (icfull, precis, icdig));
  );"""

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
REP_ICABEL = """icabel(zz) = {
  local(h, A, lv);
  if (swon && swidx>0,
    /* exp-048: the sample's radius level is a property of the cached walk
       endpoint, so derive it once per grid stretch. */
    lv = 0;
    if (icdlon,
      lv = swlev[swidx];
      if (lv == 0,
        lv = 1 + floor(icnlev * abs(zz-circc) / icrmax);
        if (lv < 1, lv = 1);
        if (lv > icnlev, lv = icnlev);
        swlev[swidx] = lv;
      );
    );
    if (icfull,
      A = rlnlm*(log(I*(zz-L))-Pi*I/2) + rlnlm2*(log(-I*(zz-L2))+Pi*I/2) + sfunczero;
      h = subst(if (lv, icdl[lv], ct), x, (zz-circc));
      icA[swidx] = A;
      icvals[swidx] = h;
    ,
      h = icvals[swidx] + subst(if (lv, icdl[lv], icdct), x, precision(zz-circc, icdig));
      icvals[swidx] = h;
      A = icA[swidx];
    );"""

PATCHES = [
    (PAT_GLOBALS, REP_GLOBALS),
    (PAT_FUNC, REP_FUNC),
    (PAT_ALLOC, REP_ALLOC),
    (PAT_BUILD, REP_BUILD),
    (PAT_ICABEL, REP_ICABEL),
]


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
    print(f"wrote {build(a.levels, HERE / f'_exp048_{a.levels}.gp')} (levels={a.levels})")
