"""Build a GMRES-accelerated engine (exp-070a) from the CURRENT fork.

Writes research/tools/_gmres.gp -- a TEST artifact, not the shipped engine.
The shipped fork stays untouched until this is validated against the gate.

Why GMRES and not Anderson: same mathematics on an affine fixed point, but the
Anderson prototype solves its least squares through the normal equations
DF^T DF, and those go singular at m ~ 19 ("impossible inverse in gauss").
Arnoldi orthogonalises explicitly and ran to k = 34 at the same state.

WHY THE HOOK SITS WHERE IT DOES -- two failed attempts, both measured:

  1. Hooking after the loop's own renormalisation (before staylor) makes andG
     renormalise twice. Result: 80x faster, sexp(0.5) wrong from the 9th digit,
     and the roundtrip still read 1.4e-211. K0 (roundtrip) is NECESSARY BUT NOT
     SUFFICIENT: sexp(slog(z)) = A^-1(A(z)) = z for ANY self-consistent Abel
     function, whatever theta it belongs to. Only agreement with the frozen
     reference is a correctness signal.

  2. Reordering the operator to y-coordinates instead of moving the hook gives
     Ghat = theta . renorm . staylor, which is NOT AFFINE: measured
     |A(2v)-2A(v)|/|A(v)| = 0.63 against 1.03e-168 for the original G. The theta
     update is not a map on ct at all -- it is a SIDE EFFECT that sets the global
     series `tht`, which staylor then reads. So the order matters, and only
     staylor . theta . renorm is affine.

Hence: the hook goes at the TOP of the loop body, and the operator keeps the
original order. Because a burst only runs while the grid is stable, the
previous iteration's `ctsamples`/`thsamples` are the right values to use, which
also freezes thsamples for the burst exactly as affstep freezes affths.

`andon` defaults to 0, so a patched engine is bit-identical to the fork until
the switch is thrown. `andaffchk()` returns the affinity residual and must be
run in the same session as any timing claim.

VORBEHALT (2026-07-30): die hier erzeugte Testkopie enthaelt WEDER exp-071
(safefs/sub-eta) NOCH exp-074 (Abschneide-Schwelle). Sie ist als A/B in sich
geschlossen und die damit gemessenen Zahlen bleiben gueltig, aber sie ist
nicht der ausgelieferte Fork. Wer sie fuer neue Messungen benutzt, misst einen
Stand von vor 0.1.1.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = Path(__file__).resolve().parent / "_gmres.gp"

GLOBALS = """quietmode=0;
andon=0; andmax=120; andcap=0.1998; andr=0; andcts=0; andths=0;
andstop=0; andlast=0; andtarget=0; andsteps=0; andused=0; andskip=0; andcts0=0;
andburst=0; andtrace=0; andgeo=0; andgmax=0; andq=0; andprec=0; andrelax=0; andgz=0; andzero=0;
/* I added || (real(Period)>47)"""

OPERATOR = r"""/* exp-070a: one Picard step as an operator on a coefficient vector, at the
   grid the loop is currently using. SAME ORDER as affstep (renorm, theta,
   staylor) -- any other order is not affine, see the module docstring. */
andtov(p, D) = vector(D, i, polcoeff(p, i-1));
/* native PARI vector ops: the element-wise GP loops cost 2.4x per operator
   application at D ~ 1300 -- the operation count was never the problem, the
   interpreter was. */
anddot(a, b) = a * b~;
andnrm(a) = sqrt(a * a~);

andG(v) = {
  local(z, sav);
  sav = default(realprecision);
  if (andrelax && (andprec > 30) && (andprec < sav), default(realprecision, andprec));
  ct = Polrev(v);
  ct = precision(ct, if (andrelax && (andprec > 30), min(precis, andprec), precis));
  ct = ct + renormr(ct);
  if (thetamode,
    thetamode++;
    tht = thtaylor(1, andths);
    if (complextaylor, tht2 = thtaylor(2, andths));
  );
  z = staylor(circc, andr, andcts);
  ct = z[1];
  andstop = z[2];
  default(realprecision, sav);
  return(andtov(ct, #v));
}

/* v - (G(x0+v) - g0), i.e. the action of (I - A) */
andop(v, g0, x0) = {
  local(a);
  if (andzero, a = andG(v); return(v - a + andgz));
  a = andG(x0 + v);
  return(v - a + g0);
}

/* MANDATORY control: GMRES is meaningless unless A is affine. Deterministic
   alternating perturbation so it is not merely a multiple of x0. */
andaffchk() = {
  local(D, x0, g0, v1, r1, r2, a1, a2);
  D = andcts + 8;
  x0 = andtov(ct, D);
  g0 = andG(x0);
  v1 = vector(D, i, x0[i]*1e-30*(-1)^i);
  r1 = andG(x0 + v1);
  r2 = andG(x0 + 2*v1);
  a1 = r1 - g0;
  a2 = r2 - g0;
  if (andnrm(a1) == 0, return(-1));
  return(andnrm(a2 - 2*a1)/andnrm(a1));
}

/* GMRES on (I-A)d = g0 - x0 at the current grid, stopped at GRID CAPACITY
   (andcts*andcap digits), not at the residual floor: past capacity it converges
   accurately to the solution of the too-coarsely discretised problem. */
gmsolve() = {
  local(D, x0, g0, bv, bet, V, H, kk, i, j, w, tol, Hk, Hs, rhs, y, yn, xn, rn, rv, cc);
  D = andcts + 8;
  andburst = 1;
  andprec = precis;
  x0 = andtov(ct, D);
  g0 = andG(x0);
  /* G(0) only AFTER G(x0): staylor caches the walk endpoints, and a walk taken
     from the zero polynomial locks in different branches -- which destroys the
     very branch-locking that makes the operator affine. Measured with the wrong
     order: 8 of 400 digits and the grid oscillating at 64. */
  if (andzero, andgz = andG(vector(D, i, 0)));
  bv = g0 - x0;
  bet = andnrm(bv);
  if (bet == 0, ct = Polrev(g0); andused = 1; andburst = 0; return(1));
  tol = bet * 10^(-andtarget);
  V = vector(andmax+1);
  H = matrix(andmax+1, andmax);
  V[1] = bv/bet;
  kk = 0; y = 0;
  for (j=1, andmax,
    w = andop(V[j], g0, x0);
    for (i=1, j,
      H[i,j] = anddot(V[i], w);
      w = w - H[i,j]*V[i];
    );
    for (i=1, j,
      cc = anddot(V[i], w);
      H[i,j] = H[i,j] + cc;
      w = w - cc*V[i];
    );
    H[j+1,j] = andnrm(w);
    /* lucky breakdown: the Krylov space is exhausted. Solving the (j+1)xj
       least squares through the normal equations then hits a singular matrix
       ("impossible inverse in gauss") -- the check has to come BEFORE the
       solve, not after. Solve the square system instead and stop. */
    if (H[j+1,j] <= bet*1e-200,
      Hs = matrix(j, j, a, b, H[a,b]);
      rhs = vector(j, a, if (a==1, bet, 0))~;
      yn = iferr(matsolve(Hs, rhs), E, 0);
      if (type(yn) == "t_COL", y = yn; kk = j);
      break;
    );
    Hk = matrix(j+1, j, a, b, H[a,b]);
    rhs = vector(j+1, a, if (a==1, bet, 0))~;
    yn = iferr(matsolve(Hk~*Hk, Hk~*rhs), E, 0);
    if (type(yn) != "t_COL", break);
    y = yn; kk = j;
    rv = rhs - Hk*y;
    rn = sqrt(sum(a=1, j+1, rv[a]^2));
    if (rn < tol, break);
    andprec = ceil(andtarget + log(rn/bet)/log(10)) + 40;
    if (andprec > precis, andprec = precis);
    if (andprec < 30, andprec = 30);
    V[j+1] = w/H[j+1,j];
  );
  if (kk == 0, andG(x0); andused = 0; andburst = 0; return(0));
  xn = x0;
  for (i=1, kk, xn = xn + y[i]*V[i]);
  /* one application lands on the image (what the loop expects next) and yields
     a consistent stopterms. It MUST run at full precision: inexact Krylov
     licenses sloppy MATVECS (they only steer the search direction), never the
     final solution. Relaxing it too gave 8x and 58 of 400 digits, and the run
     stalled at grid 512 because re could not advance. */
  andprec = precis;
  andG(xn);
  andburst = 0;
  andsteps = andsteps + kk + 2;
  andused = kk;
  if (andtrace, print("BURST cts=", andcts, " target=", precision(andtarget,6), " k=", kk, " drop=", precision(-log(rn/bet)/log(10), 6)));
  return(kk);
}

sexpinit(b,nlim,nskip,looplim) = {"""

HOOK_OLD = """    ct=precision(ct,precis);
    ct=ct+rr;
    if (thetamode,
      thetamode++;  /* used by renormr */
      /* exp-063: tht_re_mult starts NEGATIVE (log(10)/log(7/2000)) and only
         turns positive once a nonzero theta coefficient is seen. If that
         coefficient were ever exactly zero the seed would survive, and past
         re ~ 14.7 this expression goes negative straight into vector(). 6 is
         a length the code already uses on purpose (thtaylor(1,6)). */
      thsamples=max(6, floor(re*tht_re_mult+6));
      tht=thtaylor(1,thsamples);
      if (complextaylor, tht2=thtaylor(2,thsamples));  /* might need thsamples2 ... */
      if (n<=2, m=1, m=floor(thsamples/2));
      thlog=abs(polcoeff(tht,m));
      if (thlog<>0,
        thlog=log(thlog)+0.69;           /* 2x = +0.69 buffer for tht_re_mult */
        tht_re_mult=-m*1.05*log10/thlog; /* this is the multiplier for re; 1.05x is a buffer */
      );
    );

    ctsamples = floor((stopterms+20)*1.03);
    ctsamples = 4*floor(ctsamples/4);
    if ((re<relast) || (re>=skipdec),
      nskip--;
      if (re<skipdec,
        if (ctsamples<length(ct), ctsamples=length(ct)-1+4);
      );
    );
    if (ctsamples<length(ct), ctsamples=length(ct)-1);
    z=staylor(circc,r,ctsamples); ct=z[1]; stopterms=z[2];"""

# thsamples and ctsamples are hoisted ABOVE the branch. In the second attempt
# they sat inside the guarded block, so after a burst ctsamples was never
# recomputed, andlast stayed equal, and every later iteration re-solved the
# same already-solved system: re froze at 79.32490186283808127 for six
# iterations and the grid never grew (80 correct digits instead of 200).
HOOK_NEW = """    if (thetamode, thsamples=max(6, floor(re*tht_re_mult+6)));
    ctsamples = floor((stopterms+20)*1.03);
    ctsamples = 4*floor(ctsamples/4);
    if ((re<relast) || (re>=skipdec),
      nskip--;
      if (re<skipdec,
        if (ctsamples<length(ct), ctsamples=length(ct)-1+4);
      );
    );
    if (ctsamples<length(ct), ctsamples=length(ct)-1);
    if (andon && andgeo && (andlast > 0),
      /* Nested iteration: hold the level until its capacity is spent, THEN
         double. Doubling on every growth request instead ran 20 -> 52 -> 64 ->
         128 -> 256 -> 640 -> 1044 -> 1280 in eight iterations while re was
         still 14, so the first burst started COLD on the most expensive grid.
         The quadratic curve was measured from a WARM state; from far away the
         residual has broad spectral content and GMRES converges slowly
         (measured: k=61 for 35 digits, against k=17 for 63 once warm). */
      if (re >= andlast*andcap*0.9,
        ctsamples = 2*andlast;
        andgmax = 4*ceil(looplim/andcap/4) + 8;
        if (ctsamples > andgmax, ctsamples = andgmax);
      ,
        ctsamples = andlast;
      );
      if (ctsamples < length(ct), ctsamples = length(ct)-1);
    );
    if (andon && andgeo,
      andq = ctsamples*andcap;
      if (andq > precis, andq = precis);
      if (andq > looplim, andq = looplim);
      default(realprecision, max(48, min(precis, 2*floor(andq) + 60)));
    );
    andskip = 0;
    if (andon && (n > 2) && (andlast == ctsamples), andskip = 1);
    andlast = ctsamples;
    if (andskip,
      /* The theta grid is the OTHER binding resource. Measured with the
         coupled ladders: at thsamples 11 a burst got 42 digits out of 58
         steps; at thsamples 27 it got 128 out of 26 -- the quadratic regime.
         An earlier attempt sized it via tht_re_mult and failed because that
         multiplier is not yet calibrated early in the run (it is fitted from
         the theta coefficients each iteration). The trace gives the relation
         empirically instead: thsamples ~ 0.35*re, consistent over 27/79,
         39/118, 55/156, 71/205. Size it for the TARGET. */
      andr = r; andcts = ctsamples; andths = thsamples;
      andtarget = andcts*andcap;
      if (andtarget > precis, andtarget = precis);
      if (andtarget > looplim, andtarget = looplim);
      andtarget = andtarget - re;
      if (andtarget < 1, andtarget = 1);
      if (thetamode && andgeo, andths = max(andths, ceil(0.35*(re + andtarget)) + 6));
      gmsolve();
      stopterms = andstop;
    ,
      ct=precision(ct,precis);
      ct=ct+rr;
      if (thetamode,
        thetamode++;
        tht=thtaylor(1,thsamples);
        if (complextaylor, tht2=thtaylor(2,thsamples));
        if (n<=2, m=1, m=floor(thsamples/2));
        thlog=abs(polcoeff(tht,m));
        if (thlog<>0,
          thlog=log(thlog)+0.69;
          tht_re_mult=-m*1.05*log10/thlog;
        );
      );
      z=staylor(circc,r,ctsamples); ct=z[1]; stopterms=z[2];
    );"""

ICFULL_OLD = '    if ((ickey == [samples, w, r]) && (type(icct) == "t_POL"),'
ICFULL_NEW = '    if ((andburst == 0) && (ickey == [samples, w, r]) && (type(icct) == "t_POL"),'

PATCHES: list[tuple[str, str]] = [
    ("quietmode=0;\n/* I added || (real(Period)>47)", GLOBALS),
    (HOOK_OLD, HOOK_NEW),
    (ICFULL_OLD, ICFULL_NEW),
    ("sexpinit(b,nlim,nskip,looplim) = {", OPERATOR),
]


def main() -> int:
    text = SRC.read_text(encoding="utf-8")
    for i, (pat, rep) in enumerate(PATCHES):
        n = text.count(pat)
        if n != 1:
            print(f"patch {i} matched {n} times -- abort", file=sys.stderr)
            return 1
        text = text.replace(pat, rep, 1)
    DST.write_bytes(text.encode("utf-8"))
    print(f"wrote {DST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
