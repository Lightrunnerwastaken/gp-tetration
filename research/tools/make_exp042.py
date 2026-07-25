"""exp-042 candidates: fix thtaylor's extraction block.

Section profile (research/tools/section_profile.py, base e): thtaylor is
29.5% of sexpinit at dps 200 and 16.7% at dps 400 -- a block that neither
cost model accounts for. Its extraction is still what staylor lost in
exp-011b/020:

  * a quadratic rotation DFT (samples x terms complex mults at FULL precision)
  * `wtaylor = wtaylor + tot*x^s` inside the loop, i.e. O(n^2) coefficient
    copies through n freshly allocated polynomials

Two variants so the effects are separable:

  a  Polrev only     -- collect coefficients in a vector, build the polynomial
                        once. Bit-identical arithmetic; pure allocation win.
  b  Polrev + FFT    -- also replace the rotation DFT by the same transform
                        staylor uses:
                          coeff_s = (1/n) * sum_t t_est[t] * conj(tcrc[t])^s
                                  = (1/n) * (conj(c0)/om)^s * G[(s%n)+1]
                        with om = exp(2*Pi*I/n), c0 = exp(-Pi*I*(1+1/n)),
                        G = fft(powers(om^-1, n-1), t_est)  (bluedft when n
                        is not a power of two).

exp-025 tried the DFT half alone on 2026-07-11 and measured it neutral -- but
that was before exp-026/027 changed thtaylor's cost balance, so it is retested
here together with the allocation fix.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
HERE = Path(__file__).resolve().parent

PAT_LOCAL = ("thtaylor(n,samples) = {\n"
             "  local(s,t,x1,y,z,tot,t_est,tcrc,halfsamples,wtaylor,terms,thsc);")
REP_LOCAL = ("thtaylor(n,samples) = {\n"
             "  local(s,t,x1,y,z,tot,t_est,tcrc,halfsamples,wtaylor,terms,thsc,"
             "om,c0,G,cf,pw);")

PAT_EXTRACT = """  for (s=0,terms,
    tot=0;
    for (t=1,samples,
      tot=tot+t_est[t];
      t_est[t]=t_est[t]*conj(tcrc[t]);
    );
    tot=tot/samples;
/*  if (s>=1, tot=tot*(rinv)^s );   */
    wtaylor=wtaylor+tot*x^s;
  );
  wtaylor=precision(wtaylor,precis);"""

REP_A = """  /* exp-042a: same arithmetic, but collect the coefficients and build the
     polynomial once -- the old form allocated a fresh polynomial per term
     (O(n^2) coefficient copies at full precision). */
  cf = vector(terms+1, i, 0);
  for (s=0,terms,
    tot=0;
    for (t=1,samples,
      tot=tot+t_est[t];
      t_est[t]=t_est[t]*conj(tcrc[t]);
    );
    cf[s+1]=tot/samples;
  );
  wtaylor=Polrev(cf);
  wtaylor=precision(wtaylor,precis);"""

REP_B = """  /* exp-042b: replace the quadratic rotation DFT with the transform
     staylor has used since exp-011b/020, and build the polynomial once.
       coeff_s = (1/n) sum_t t_est[t] conj(tcrc[t])^s
               = (1/n) (conj(c0)/om)^s G[(s%n)+1]
     om = exp(2*Pi*I/n), c0 = exp(-Pi*I*(1+1/n)),
     G = fft(powers(om^-1, n-1), t_est); bluedft for non-power-of-two n. */
  om = exp(2*Pi*I/samples);
  c0 = exp(-Pi*I*(1+1/samples));
  if (2^valuation(samples,2)==samples,
    G = fft(powers(om^(-1), samples-1), t_est)
  ,
    G = bluedft(t_est));
  pw = powers(conj(c0)/om, terms);
  cf = vector(terms+1, j, pw[j] * G[((j-1)%samples)+1] / samples);
  wtaylor=Polrev(cf);
  wtaylor=precision(wtaylor,precis);"""


def build(variant: str, dst: Path) -> Path:
    text = SRC.read_text(encoding="utf-8")
    rep = {"a": REP_A, "b": REP_B}[variant]
    for pat, r in ((PAT_LOCAL, REP_LOCAL), (PAT_EXTRACT, rep)):
        if text.count(pat) != 1:
            print(f"pattern matched {text.count(pat)} times -- abort", file=sys.stderr)
            raise SystemExit(1)
        text = text.replace(pat, r, 1)
    dst.write_bytes(text.encode("utf-8"))
    return dst


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--variant", choices=["a", "b"], default="b")
    a = ap.parse_args()
    p = build(a.variant, HERE / f"_exp042{a.variant}.gp")
    print(f"wrote {p}")
