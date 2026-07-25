"""exp-056 candidate: exact-N DFT by radix-r decimation instead of Bluestein.

The extraction is 11.3% of the run at dps 400 and routes every non-power-of-two
grid through bluedft(), which needs FFTs of length >= 3N -- 4 to 8 times the
data length. But the grids are never arbitrary: the quantum is 256 = 2^8, so
`samples` is always m*2^8 with m = 1..9, i.e. of the form r*2^k with an odd
r <= 9:

    samples  factorization   bluedft FFT length   radix split
      768     3 * 2^8              8192           3 FFTs of 256
     1280     5 * 2^8              8192           5 FFTs of 256
     1792     7 * 2^8             16384           7 FFTs of 256
     2304     9 * 2^8             16384           9 FFTs of 256

Decimation in time: write j-1 = r*a + s with a in [0,m), s in [0,r), and use
omega^r = exp(-2*Pi*I/m):

    X[k] = sum_s omega^(s*k) * SUB_s[k mod m],
    SUB_s = DFT_m( t[s+1], t[r+s+1], t[2r+s+1], ... )

so r power-of-two FFTs of length m = N/r plus N*r twiddle multiplications. For
the real-side call (which doubles the vector to 4608 = 9*512) that is roughly
82k complex multiplications against Bluestein's ~688k -- about 8x.

The split is an exact reorganization of the same sum, so it cannot be less
accurate than Bluestein; it avoids the chirp vectors entirely. Anything that
does not factor as r*2^k with r <= 9 still falls through to bluedft().
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = Path(__file__).resolve().parent / "_exp056.gp"

PAT_FUNC = "bluedft(t) = {"
REP_FUNC = """/* exp-056: exact-N DFT for N = r*2^k with small odd r, by radix-r decimation
   into r power-of-two FFTs plus twiddles. Same convention as bluedft:
   X[k] = sum_j t[j] * omega^((j-1)(k-1)),  omega = exp(-2*Pi*I/N).
   Falls back to bluedft for anything that does not factor that way. */
mixdft(t) = {
  local(nn, m, r, w, sub, s, a, j, omp, res, idx, acc);
  nn = #t;
  m = 2^valuation(nn, 2);
  r = nn/m;
  if ((r > 9) || (m < 8), return(bluedft(t)));
  w = powers(exp(-2*Pi*I/m), m-1);
  sub = vector(r, s, fft(w, vector(m, a, t[(a-1)*r + s])));
  if (r == 1, return(sub[1]));
  omp = powers(exp(-2*Pi*I/nn), nn-1);
  res = vector(nn);
  for (j=1, nn,
    idx = ((j-1) % m) + 1;
    acc = sub[1][idx];
    for (s=2, r,
      acc = acc + omp[(((j-1)*(s-1)) % nn) + 1] * sub[s][idx];
    );
    res[j] = acc;
  );
  return(res);
}

bluedft(t) = {"""

CALLS = [
    ("    G = bluedft(t_est));", "    G = mixdft(t_est));"),
    ("        G = bluedft(exin));", "        G = mixdft(exin));"),
    ("        G = bluedft(concat(exin, vector(samples, i, 0))));",
     "        G = mixdft(concat(exin, vector(samples, i, 0))));"),
]


def main() -> int:
    text = SRC.read_text(encoding="utf-8")
    if text.count(PAT_FUNC) != 1:
        print(f"function anchor matched {text.count(PAT_FUNC)} times -- abort", file=sys.stderr)
        return 1
    text = text.replace(PAT_FUNC, REP_FUNC, 1)
    for pat, rep in CALLS:
        if text.count(pat) != 1:
            print(f"call site {pat!r} matched {text.count(pat)} times -- abort", file=sys.stderr)
            return 1
        text = text.replace(pat, rep, 1)
    DST.write_bytes(text.encode("utf-8"))
    print(f"wrote {DST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
