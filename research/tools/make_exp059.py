"""exp-059 candidate: the sampling grid is a geometric sequence, not N exps.

Both sampling loops build their grid point by point with a full-precision
complex exponential:

  staylor, complex path :  x1 = -1 + -1/samples + 2*s/samples ; tcrc[s]=exp(Pi*I*x1)
  staylor, real path    :  x1 = -1/(2*samples) + s/samples    ; tcrc[s]=exp(Pi*I*x1)
  thtaylor              :  x1 = -1 + -1/samples + 2*s/samples ; tcrc[s]=exp(Pi*I*x1)

x1 is affine in s, so tcrc[s] = c * mu^s with mu and c fixed -- and the very
same mu and c are rebuilt 20 lines below for the extraction. The loop spends N
complex exponentials on a geometric sequence.

A naive `powers(mu, N)` would accumulate ~N ulp (reproduced elsewhere: ~2000 ulp
at N=1024/dps 404, which eats 3 digits of headroom). This uses a DOUBLING build
instead -- log2(N) exponentials and a table combination -- so the error stays
O(log N) ulp while the cost stays O(N) multiplications.

Stacked with it: thtaylor's `y = log(z)/(2*Pi*I)` in thfunc undoes the exp it
just did. For x1 strictly inside (-1,1) -- true for every s -- the principal log
gives exactly Pi*I*x1, so y = x1/2, and x1 is an exact t_FRAC, making the direct
value strictly MORE accurate than the round trip. (Tried alone as exp-050 and
measured below the noise floor; it belongs with this change, not against it.)

Neither part changes the mathematics, so the anchor must stay bit-identical.
That, not the clock, is the acceptance test -- both effects are a few percent
and the frozen benchmark's noise threshold is 6.8%.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = Path(__file__).resolve().parent / "_exp059.gp"

PAT_HELPER = "staylor( w,r,samples) = {"
REP_HELPER = """/* exp-059: tcrc[s] = c*mu^s built by doubling. `powers(mu,n)` would carry
   ~n ulp; a doubling table carries ~log2(n) ulp for the same O(n) work. */
geoseq(c, mu, n) = {
  local(v, half, j, blk);
  v = vector(n);
  if (n <= 0, return(v));
  v[1] = c*mu;
  blk = 1;
  while (blk < n,
    half = min(blk, n-blk);
    for (j=1, half, v[blk+j] = v[blk] * v[j] / c);
    blk = blk + half;
  );
  return(v);
}

staylor( w,r,samples) = {"""

PAT_CPLX = """    for(s=1, samples,
      x1=-1+-1/(samples)+(2*s/samples);
      tcrc[s]=exp(Pi*I*x1); /* -Pi to Pi */
      swidx = s;
      t_est[s]=sfunc(w+r*tcrc[s]*exp(I*argc))[1];
    );"""
REP_CPLX = """    /* exp-059: x1 is affine in s, so tcrc is geometric */
    tcrc = geoseq(exp(Pi*I*(-1-1/samples)), exp(2*Pi*I/samples), samples);
    y0 = exp(I*argc);
    for(s=1, samples,
      swidx = s;
      t_est[s]=sfunc(w+r*tcrc[s]*y0)[1];
    );"""

PAT_REAL = """    for(s=1, samples,
      x1=-1/(2*samples)+(s/samples);
      tcrc[s]=exp(Pi*I*x1);
      swidx = s;
      t_est[s]=sfunc(w+r*tcrc[s])[1];
    );"""
REP_REAL = """    /* exp-059: x1 is affine in s, so tcrc is geometric */
    tcrc = geoseq(exp(-Pi*I/(2*samples)), exp(Pi*I/samples), samples);
    for(s=1, samples,
      swidx = s;
      t_est[s]=sfunc(w+r*tcrc[s])[1];
    );"""

PAT_TH = """  for(s=1, samples,
    x1 = -1 + -1/(samples) + (2*s/samples); /* -Pi to Pi */
    tcrc[s] = exp(Pi*I*x1);
    thidx = s;
    t_est[s]= thfunc(tcrc[s],n);
  );"""
REP_TH = """  /* exp-059: geometric grid, and hand thfunc the exact y = x1/2 instead of
     making it recover x1 from exp(Pi*I*x1) with a full-precision log. */
  tcrc = geoseq(exp(Pi*I*(-1-1/samples)), exp(2*Pi*I/samples), samples);
  for(s=1, samples,
    x1 = -1 + -1/(samples) + (2*s/samples);
    thidx = s;
    t_est[s]= thfunc(tcrc[s], n, x1/2, 1);
  );"""

PAT_SIG = """thfunc(z,n) = {
  local(y,p,h,A);
  if (n<>2,
    y=log(z)/(2*Pi*I);"""
REP_SIG = """thfunc(z,n,yin,yon) = {
  local(y,p,h,A);
  if (n<>2,
    /* exp-059: the caller knows y exactly; this log only undoes its exp */
    y = if (yon, yin, log(z)/(2*Pi*I));"""

PAT_N2 = """    y=log(z)/(-2*Pi*I);
    y = abelest(superf2(ztl+y),ct)-y-ztl;"""
REP_N2 = """    y = if (yon, -yin, log(z)/(-2*Pi*I));
    y = abelest(superf2(ztl+y),ct)-y-ztl;"""

PATCHES = [(PAT_HELPER, REP_HELPER), (PAT_CPLX, REP_CPLX), (PAT_REAL, REP_REAL),
           (PAT_SIG, REP_SIG), (PAT_N2, REP_N2), (PAT_TH, REP_TH)]


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
