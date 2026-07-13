# gp-tetration — fast Kneser tetration with verified references and a base atlas

A Python wrapper and optimized fork of Sheldon Levenstein's (Sheldonison's)
`fatou.gp` — the classic PARI/GP implementation of Kneser's real-analytic
tetration. Everything here builds on that engine; full credit for the
underlying construction and original code goes to him.

## What this adds on top of fatou.gp

**1. An optimized engine fork (~13× faster at high precision).**
21 gate-verified optimizations (`src/fatou_backend/vendor/fatou_fork.gp`):
persistent walk/Schröder caches, incremental Taylor/theta/extraction updates
at reduced precision, FFT/Bluestein extraction on exact-size grids, precision
laddering, contour-radius tuning. Every change had to pass a frozen accuracy
gate (all bases, full digits) before being kept; methods, measurements and
the negative results (what did *not* work, and why) are documented in
`research/METHODS.md`.

| target (base e) | true digits | original | fork |
|---|---|---|---|
| dps 300 | 294 | ~27 min | **~3 min** |
| dps 520 | 495 | ~5.4 h | **~25 min** |
| dps 1020 | 973 | (memory crash) | **~6.7 h** |

**2. Reference values with *proven* error bounds.**
`research/reference/values.json` contains sexp values whose accuracy is
established by the error-vector method (two independent runs at different
iteration depths; their agreement bounds the worse run's true error) and,
for the deep tiers, additionally by engine diversity (fork vs. original):

- `sexp_e(0.5)` to **972 proven digits** (dps 1020/1033 pair)
- `sexp_e(0.5)` to 698 / 497 proven digits (700 / 500 tiers)
- `sexp_2(0.5)` to 497 proven digits (error vector **and** engine diversity)

**3. A base atlas: tetration for any base from one anchor.**
`fatou_backend.basechange` implements the base-change ladder: after a single
base-e setup, `sexp_anchor` / `slog_anchor` evaluate tetration for *any* base
via a small Fourier-mode table (`research/reference/phi_modes.json`,
on-demand, ~3 KB per base, first request ≈ 0.2 s):

- validated against the proven references: ~1e-24 absolute (k=20 head,
  dps-60 table); ~1e-37 with the k=40/dps-113 table; scales further
- round-trips `slog(sexp(y))` consistent to ~1e-52
- certified-enclosure prototype (ball arithmetic) in
  `research/tools/m54_cert_proto.py`

The mathematics (base-change phase law, mode decay, universality, µ-hub)
is developed in the mixed-base tetration paper series — see `papers/` —
and this repository is its computational companion.

**4. A zero-dependency interactive demo.**
`research/tools/tetration_calculator.html` — a single static HTML file that
plots sexp_b(x) with a *continuous* base slider (1.5 … 100), powered by the
atlas tables in plain JavaScript. No server, no PARI. Typical accuracy
1e-6…1e-9 (less near the base edge).

## Install / use

Requirements: PARI/GP ≥ 2.15 (`gp_exe=` parameter or `FATOU_GP_EXE`),
Python ≥ 3.11, `mpmath` (optionally `python-flint` for certification tools).

```python
from fatou_backend.gp_backend import FatouGP
from fatou_backend import basechange
import mpmath as mp

gp = FatouGP(dps=60)                      # persistent worker + state cache
gp.sexp(2, mp.mpf("0.5"))                 # direct engine call (base 2)
basechange.sexp_anchor(gp, 3, "0.5")     # any base via the e-anchor, ~ms
```

- `fatou.gp` resolution: `fatou_gp=` parameter → `FATOU_GP_FILE` →
  vendored file. The optimized fork: `fatou_gp=".../vendor/fatou_fork.gp"`.
- First initialization of a (base, dps) pair is computed once and cached in
  `~/.cache/fatou_backend/`; later sessions start in ~0.1 s.
- Batch/CLI: `python -m fatou_backend.cli sexp --base e --values 0.5`;
  `FatouGP(n_workers=N)` parallelizes large batches.

Tests: `python -m pytest tests/` (53 tests, includes atlas validation
against the proven references; set `FATOU_BACKEND_RUN_SLOW=1` for the
slow gate block).

## Reproducibility

- `bench/` — frozen workload/metrics benchmark (`bench/results/`)
- `research/gate.py` + `research/reference/` — the frozen accuracy gate used
  for every keep (anti-gaming: gate and references are immutable)
- `research/METHODS.md` — verification method, the full keep stack,
  measured speedups, and all negative results
- `papers/` — the mixed-base tetration paper series (the mathematics
  behind the base atlas)

## Credits & license

- `fatou.gp` © Sheldon Levenstein, published on the Tetration Forum;
  vendored unmodified as `src/fatou_backend/vendor/fatou.gp`, optimized fork
  alongside as `fatou_fork.gp`. This project exists thanks to that work.
- Base-change mathematics: the mixed-base tetration paper series
  (J. Justus, 2026), included in `papers/`.
- Wrapper, fork optimizations, references, atlas: this repository, 2026.
