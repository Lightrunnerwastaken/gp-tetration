"""exp-049 candidate: the same per-radius truncation for thtaylor's Horner block.

exp-048 truncates ct per sample radius in staylor's sampling loop. thtaylor
has its own, separate Horner block (thfunc -> subst(thdct, ...)) which the
section profile puts at 21.4% of the run at dps 400, and it still evaluates
the full degree.

Its points are cheaper to handle than staylor's, because they sit in a NARROW
band. Measured at dps 300 (base e):

    |p - circc| over the 192 theta samples:
        min 0.567   q1 0.579   median 0.616   q3 0.673   max 0.747
    circr = 1.337, thdig = 159, deg(thdct) = 1792

    terms needed at the OUTER edge: 159 / log10(1.337/0.747) = 629
    -> one truncation level covers every theta sample, ~2.8x fewer terms

So this needs no bucketing at all: track the largest radius seen on the grid
(the points are cached in thsw and fixed per stretch), and truncate once per
pass at that radius plus a 2% margin.

Safety: the first pass of a stretch runs with thfull=1 and the full polynomial,
which is also when the radius gets measured -- so the truncation is only ever
built from a radius that has actually been observed on the current grid.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = Path(__file__).resolve().parent / "_exp049.gp"

PAT_GLOBALS = "icnlev=8; swlev=0; icdl=0; icrmax=1; icdlon=0;"
REP_GLOBALS = ("icnlev=8; swlev=0; icdl=0; icrmax=1; icdlon=0;\n"
               "/* exp-049: same idea on the theta side -- one level suffices there */\n"
               "thdl1=0; thrmax=0; thdlon=0;")

# single-radius truncation helper, placed next to icbuild
PAT_HELPER = "icbuild(pp, dig) = {"
REP_HELPER = """/* exp-049: truncate pp to the degree needed at radius rho, same criterion
   as icbuild but for a single radius (the theta points sit in a narrow band). */
ictrunc(pp, dig, rho) = {
  local(v, n, lg2, tol2, lr2, kk, ii, mx);
  if (type(pp) != "t_POL", return(pp));
  v = Vecrev(pp);
  n = #v;
  if (n < 64, return(pp));
  lg2 = vector(n, kk, if (v[kk]==0, -1000000, exponent(v[kk])));
  mx = vecmax(lg2);
  if (mx <= -1000000, return(pp));
  tol2 = mx - dig*3.3219280948873623 - 16;
  lr2 = log(rho)/log(2);
  kk = n;
  while ((kk > 1) && (lg2[kk] + (kk-1)*lr2 < tol2), kk--);
  return(if (kk >= n, pp, Polrev(vector(kk, ii, v[ii]))));
}

icbuild(pp, dig) = {"""

PAT_ALLOC = """    if (thfull,
      thsw = vector(samples); thswv = vector(samples);
      thicv = vector(samples); thicA = vector(samples);
    );
    thskey = samples;
    thon = 1;
  );"""
REP_ALLOC = """    if (thfull,
      thsw = vector(samples); thswv = vector(samples);
      thicv = vector(samples); thicA = vector(samples);
      thrmax = 0;   /* exp-049: radius not yet observed on this grid */
    );
    thskey = samples;
    thon = 1;
    /* exp-049: the theta samples sit in a narrow radius band, so one
       truncation of the diff polynomial covers all of them. Only built once
       the grid's radius has actually been observed (i.e. from the second pass
       of a stretch on); the first pass runs full and measures it. */
    thdlon = 0;
    if ((thfull == 0) && (thrmax > 0),
      thdl1 = ictrunc(thdct, thdig, thrmax*1.02);
      if (poldegree(thdl1) < poldegree(thdct), thdlon = 1);
    );
  );"""

PAT_THFUNC = """      if (thfull,
        A = rlnlm*(log(I*(p-L))-Pi*I/2) + rlnlm2*(log(-I*(p-L2))+Pi*I/2) + sfunczero;
        h = subst(ct, x, (p-circc));
        thicA[thidx] = A; thicv[thidx] = h;
      ,
        h = thicv[thidx] + subst(thdct, x, precision(p-circc, thdig));
        thicv[thidx] = h;
        A = thicA[thidx];
      );"""
REP_THFUNC = """      /* exp-049: observe the radius so the next pass can truncate */
      if (abs(p-circc) > thrmax, thrmax = abs(p-circc));
      if (thfull,
        A = rlnlm*(log(I*(p-L))-Pi*I/2) + rlnlm2*(log(-I*(p-L2))+Pi*I/2) + sfunczero;
        h = subst(ct, x, (p-circc));
        thicA[thidx] = A; thicv[thidx] = h;
      ,
        h = thicv[thidx] + subst(if (thdlon, thdl1, thdct), x, precision(p-circc, thdig));
        thicv[thidx] = h;
        A = thicA[thidx];
      );"""

PATCHES = [
    (PAT_GLOBALS, REP_GLOBALS),
    (PAT_HELPER, REP_HELPER),
    (PAT_ALLOC, REP_ALLOC),
    (PAT_THFUNC, REP_THFUNC),
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
