"""exp-050 candidate: thtaylor's exp/log round trip is an identity.

thtaylor builds, per sample,

    x1      = -1 + -1/samples + 2*s/samples        (exact rational)
    tcrc[s] = exp(Pi*I*x1)
    t_est[s]= thfunc(tcrc[s], n)

and the first thing thfunc does is undo it:

    y = log(z)/(2*Pi*I)

For x1 strictly inside (-1, 1) -- which it is for every s, since
x1 ranges over (-1+1/samples, 1-1/samples) -- the principal logarithm gives
log(exp(Pi*I*x1)) = Pi*I*x1 exactly, so y = x1/2. The exp and the log cancel.

That round trip costs one exp and one log per sample at FULL working precision
on every pass. Since exp-042b replaced the rotation DFT with an FFT, tcrc is
no longer used anywhere else in thtaylor either -- the vector exists purely to
be undone.

Passing x1/2 straight in is also slightly *more* accurate: x1 is an exact
rational, so y becomes exact instead of a float round trip.

n == 2 (the lower theta map) uses y = log(z)/(-2*Pi*I) = -x1/2, so the same
value serves both branches with a sign.

This was found while chasing why exp-049 (truncating thtaylor's Horner from
1792 to 885 terms, verified active) produced no speedup at all: the Horner is
not what that loop spends its time on.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = Path(__file__).resolve().parent / "_exp050.gp"

PAT_SIG = """thfunc(z,n) = {
  local(y,p,h,A);
  if (n<>2,
    y=log(z)/(2*Pi*I);"""
REP_SIG = """thfunc(z,n,yin,yon) = {
  local(y,p,h,A);
  if (n<>2,
    /* exp-050: the caller already knows y exactly (see thtaylor); the
       log() here only undoes an exp() there. */
    y = if (yon, yin, log(z)/(2*Pi*I));"""

PAT_N2 = """    y=log(z)/(-2*Pi*I);
    y = abelest(superf2(ztl+y),ct)-y-ztl;"""
REP_N2 = """    y = if (yon, -yin, log(z)/(-2*Pi*I));
    y = abelest(superf2(ztl+y),ct)-y-ztl;"""

PAT_LOOP = """  for(s=1, samples,
    x1 = -1 + -1/(samples) + (2*s/samples); /* -Pi to Pi */
    tcrc[s] = exp(Pi*I*x1);
    thidx = s;
    t_est[s]= thfunc(tcrc[s],n);
  );"""
REP_LOOP = """  for(s=1, samples,
    x1 = -1 + -1/(samples) + (2*s/samples); /* -Pi to Pi */
    thidx = s;
    /* exp-050: hand thfunc the exact y = x1/2 instead of making it recover
       x1 from exp(Pi*I*x1) with a full-precision log. Since exp-042b the
       extraction uses an FFT, so tcrc is not needed for anything else. */
    t_est[s]= thfunc(0, n, x1/2, 1);
  );"""

PATCHES = [(PAT_SIG, REP_SIG), (PAT_N2, REP_N2), (PAT_LOOP, REP_LOOP)]


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
