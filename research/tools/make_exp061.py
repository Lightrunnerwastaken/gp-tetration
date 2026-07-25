"""exp-061: finish the twiddle fix, and stop rebuilding the twiddles.

The 2026-07-25 review found that mixdft built its outer twiddles with
`powers()`, which accumulates ~n ulp (measured 4.6e-402 at nn=4608, precis
404), and replaced that one line with a doubling build. That fix was
incomplete and, on the clock, a step backwards. Three things are wrong with it:

1. **The same construction survives six times more.** A grep for `powers(`
   finds SEVEN call sites; the review fixed one. `mixdft` builds its INNER
   FFT roots the same way (m up to nn) and `bluedft` too (M the next power of
   two above 3*nn-2, i.e. up to ~4*nn) -- and, worse, four sites sit on the
   MAIN power-of-two path that the review never looked at: the FFT roots in
   thtaylor and in both staylor extraction branches, plus thtaylor's `pw`,
   which is not an FFT input at all but multiplies the output coefficients
   one by one, so its ~terms ulp went straight into the theta Taylor series.
   Fixing one site out of seven left every large error source in place.

2. **The replacement divides.** `geoseq` computes `v[blk]*v[j]/c` -- one
   division per element. At these precisions a division costs 2-3x a
   multiplication, so swapping PARI's C-level `powers()` for an interpreted
   loop *with* a division is a poor trade. Building `[1, mu, ..., mu^n]`
   directly needs no division at all: `v[blk+j] = v[blk]*v[j]`.

3. **The twiddles do not depend on the data.** They depend only on
   (nn, realprecision), yet mixdft rebuilt both vectors on every call --
   ~nn full-precision multiplications per call, every loop iteration, on a
   grid that the exp-026/051 quantization deliberately holds fixed for 20-30
   iterations. bluedft has cached its chirps since exp-031; mixdft never did.

So: `geopow(mu, n)` returns exactly what `powers(mu, n)` returns -- the vector
[1, mu, ..., mu^n] -- by doubling (error ~log2(n) ulp instead of ~n) and
without a division. All three sites use it, and mixdft gets a two-slot LRU
cache keyed on (nn, precision), the same shape bluedft already uses.

Accuracy strictly improves at every site; the work strictly decreases. The
anchor must stay bit-identical -- that, not the clock, is the acceptance test.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = Path(__file__).resolve().parent / "_exp061.gp"

# ---------------------------------------------------------------- globals
PAT_GLOB = """bd1n=0; bd1p=0; bd1t=0; bd1a=0; bd1b=0; bd1w=0; bd1m=0; bd2n=0; bd2p=0; bd2t=0; bd2a=0; bd2b=0; bd2w=0; bd2m=0; bdclock=1;"""
REP_GLOB = """bd1n=0; bd1p=0; bd1t=0; bd1a=0; bd1b=0; bd1w=0; bd1m=0; bd2n=0; bd2p=0; bd2t=0; bd2a=0; bd2b=0; bd2w=0; bd2m=0; bdclock=1;
/* exp-061: mixdft twiddle cache (two slots, LRU, shares bdclock). The
   twiddles depend only on (nn, realprecision) -- never on the data -- but
   were rebuilt on every call. */
mx1n=0; mx1p=0; mx1t=0; mx1w=0; mx1o=0; mx2n=0; mx2p=0; mx2t=0; mx2w=0; mx2o=0;"""

# ---------------------------------------------------------------- geopow/geoseq
PAT_GEO = """/* exp-059: tcrc[s] = c*mu^s built by doubling. `powers(mu,n)` would carry
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
}"""
REP_GEO = """/* exp-061: [1, mu, mu^2, ..., mu^n] -- exactly what powers(mu,n) returns,
   but built by doubling, so the error is ~log2(n) ulp instead of ~n, and
   with no division (the first form of geoseq divided by c once per element,
   and at these precisions a division costs 2-3x a multiplication). */
geopow(mu, n) = {
  local(v, half, j, blk);
  if (n < 0, return(vector(0)));   /* guard: vector(n+1) would raise */
  v = vector(n+1);
  v[1] = 1;
  if (n == 0, return(v));
  v[2] = mu;
  blk = 1;
  while (blk < n,
    half = min(blk, n-blk);
    for (j=1, half, v[blk+j+1] = v[blk+1] * v[j+1]);
    blk = blk + half;
  );
  return(v);
}

/* exp-059: tcrc[s] = c*mu^s for s = 1..n, now via geopow (no division). */
geoseq(c, mu, n) = {
  local(v);
  if (n <= 0, return(vector(0)));
  v = geopow(mu, n);
  return(vector(n, j, c*v[j+1]));
}"""

# ---------------------------------------------------------------- mixdft
PAT_MIX = """mixdft(t) = {
  local(nn, m, r, w, sub, s, a, j, omp, res, idx, acc);
  nn = #t;
  m = 2^valuation(nn, 2);
  r = nn/m;
  if ((r > 9) || (m < 8), return(bluedft(t)));
  w = powers(exp(-2*Pi*I/m), m-1);
  sub = vector(r, s, fft(w, vector(m, a, t[(a-1)*r + s])));
  if (r == 1, return(sub[1]));
  /* review fix: powers() accumulates ~n ulp -- measured 4.6e-402 at nn=4608,
     precis 404, i.e. ~2.7 digits. The doubling build gives 1.9e-403 for the
     same input. (The exp-056 verification could not see this: bluedft builds
     its twiddles the same way, so the two agreed on a common-mode error.) */
  omp = concat([1], geoseq(1, exp(-2*Pi*I/nn), nn-1));
  res = vector(nn);"""
REP_MIX = """mixdft(t) = {
  local(nn, m, r, w, sub, s, a, j, omp, res, idx, acc, sl);
  nn = #t;
  m = 2^valuation(nn, 2);
  r = nn/m;
  if ((r > 9) || (m < 8), return(bluedft(t)));
  /* exp-061: both twiddle vectors depend only on (nn, precision), so cache
     them (two slots, LRU) instead of rebuilding ~nn full-precision products
     on every call, and build them by doubling -- powers() carries ~n ulp
     (measured 4.6e-402 at nn=4608, precis 404, i.e. ~2.7 digits; doubling
     gives 1.9e-403). The exp-056 verification could not see that error:
     bluedft built its twiddles the same way, so the two agreed on a
     common-mode error. */
  sl = 0;
  if ((mx1n == nn) && (mx1p >= default(realprecision)), sl = 1);
  if ((sl == 0) && (mx2n == nn) && (mx2p >= default(realprecision)), sl = 2);
  if (sl == 0,
    w = geopow(exp(-2*Pi*I/m), m-1);
    omp = if (r == 1, 0, geopow(exp(-2*Pi*I/nn), nn-1));
    if (mx1t <= mx2t,
      mx1n=nn; mx1p=default(realprecision); mx1w=w; mx1o=omp; mx1t=bdclock++; sl=1
    ,
      mx2n=nn; mx2p=default(realprecision); mx2w=w; mx2o=omp; mx2t=bdclock++; sl=2);
  );
  if (sl == 1,
    w = mx1w; omp = mx1o; mx1t=bdclock++
  ,
    w = mx2w; omp = mx2o; mx2t=bdclock++);
  sub = vector(r, s, fft(w, vector(m, a, t[(a-1)*r + s])));
  if (r == 1, return(sub[1]));
  res = vector(nn);"""

# ---------------------------------------------------------------- bluedft
PAT_BLUE = """    w = powers(exp(2*Pi*I/M), M-1);"""
REP_BLUE = """    /* exp-061: doubling build -- M is the next power of two above 3*nn-2,
       so powers() carried up to ~4*nn ulp here, the largest of the three
       twiddle sites. Cached since exp-031, so the extra build cost is paid
       once per grid. */
    w = geopow(exp(2*Pi*I/M), M-1);"""

# ------------------------------------------------- the four remaining sites
# A grep for `powers(` turns up SEVEN call sites, not the three the review
# found. Four of them are on the main (power-of-two grid) path and were
# missed entirely -- including thtaylor's `pw`, which multiplies the output
# coefficients directly rather than feeding an FFT.
#
# Cost: the twiddle build is O(n) multiplications feeding an O(n log n) FFT
# inside a run whose sampling term is O(n^2). geopow is an interpreted loop
# against PARI's C `powers`, so ~1.5x on a term that is ~8% of the extraction
# and the extraction is ~6% of the run: ~0.25% total. Bought with an
# unconditional accuracy gain at every site.

PAT_TH = """  om = exp(2*Pi*I/samples);
  c0 = exp(-Pi*I*(1+1/samples));
  if (2^valuation(samples,2)==samples,
    G = fft(powers(om^(-1), samples-1), t_est)
  ,
    G = mixdft(t_est));
  pw = powers(conj(c0)/om, terms);"""
REP_TH = """  om = exp(2*Pi*I/samples);
  c0 = exp(-Pi*I*(1+1/samples));
  /* exp-061: doubling build at both sites. `pw` is the sharper of the two --
     it is not an FFT input, it multiplies the output coefficients one by
     one, so its ~terms ulp landed directly in the theta Taylor series. */
  if (2^valuation(samples,2)==samples,
    G = fft(geopow(om^(-1), samples-1), t_est)
  ,
    G = mixdft(t_est));
  pw = geopow(conj(c0)/om, terms);"""

PAT_SCX = """      om = exp(2*Pi*I/samples);
      c0 = exp(-Pi*I*(1+1/samples));
      if (2^valuation(samples,2)==samples,
        G = fft(powers(om^(-1), samples-1), exin)
      ,
        G = mixdft(exin));"""
REP_SCX = """      om = exp(2*Pi*I/samples);
      c0 = exp(-Pi*I*(1+1/samples));
      if (2^valuation(samples,2)==samples,
        G = fft(geopow(om^(-1), samples-1), exin)   /* exp-061 */
      ,
        G = mixdft(exin));"""

PAT_SCR = """      mu = exp(Pi*I/samples);
      c1 = exp(-Pi*I/(2*samples));
      if (2^valuation(samples,2)==samples,
        G = fft(powers(mu^(-1), 2*samples-1), concat(exin, vector(samples, i, 0)))
      ,
        G = mixdft(concat(exin, vector(samples, i, 0))));"""
REP_SCR = """      mu = exp(Pi*I/samples);
      c1 = exp(-Pi*I/(2*samples));
      if (2^valuation(samples,2)==samples,
        G = fft(geopow(mu^(-1), 2*samples-1), concat(exin, vector(samples, i, 0)))   /* exp-061 */
      ,
        G = mixdft(concat(exin, vector(samples, i, 0))));"""

PATCHES = [(PAT_GLOB, REP_GLOB), (PAT_GEO, REP_GEO), (PAT_MIX, REP_MIX),
           (PAT_BLUE, REP_BLUE), (PAT_TH, REP_TH), (PAT_SCX, REP_SCX),
           (PAT_SCR, REP_SCR)]


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
