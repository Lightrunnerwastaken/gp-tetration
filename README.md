# gp-tetration — fast Kneser tetration with verified references and a base atlas

A Python wrapper and optimized fork of Sheldon Levenstein's (Sheldonison's)
`fatou.gp` — the classic PARI/GP implementation of Kneser's real-analytic
tetration. Everything here builds on that engine; full credit for the
underlying construction and original code goes to him.

## What this adds on top of fatou.gp

**1. An optimized engine fork (~13× faster at high precision).**
27 gate-verified optimizations (`src/fatou_backend/vendor/fatou_fork.gp`):
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
Install with `pip install -e .` in the repo root (or put `src/` on
`PYTHONPATH`).

```python
from fatou_backend.gp_backend import FatouGP
from fatou_backend import basechange
import mpmath as mp

gp = FatouGP(dps=300, fatou_gp="fork")   # the optimized engine, ~13x faster
gp.sexp(2, mp.mpf("0.5"))                 # direct engine call (base 2)
basechange.sexp_anchor(gp, 3, "0.5")     # any base via the e-anchor, ~ms
```

- **Engine selection**: `fatou_gp="fork"` (optimized) or `"original"`
  (unmodified fatou.gp); a path or `FATOU_GP_FILE` also works; default is
  the vendored original.
- **Precision**: set `dps` and expect ~`dps − 24` true digits (measured
  calibration law, see `research/METHODS.md`). The defaults fully converge
  at any `dps` on both engines — a plain `FatouGP(dps=...)` call reproduces
  the published reference values; check yourself with
  `python research/tools/verify_reference.py --key "sexp|e|0.5|500"`.
- **Caching is automatic**: the first initialization of each
  (base, dps, knobs, engine) combination is computed once and stored in
  `~/.cache/fatou_backend/` (override with `FATOU_CACHE_DIR`); later
  sessions restore it in ~0.1 s. The key includes a hash of the engine
  file, so switching engines or knobs never reuses a stale state.
- **Supported bases**: real bases > e^(1/e) (Kneser construction,
  gate-verified), real bases 1 < b < e^(1/e) (regular iteration at the
  attracting fixed point, e.g. sexp_1.2(0.5) = 1.13626…, verified
  against the unmodified original engine), and complex bases (2+I,
  0.8+0.4*I, … are gate-verified; e.g. b = −1 works but has no
  independent references). Bases 0 < b < 1 are not supported.
- Batch/CLI:
  `python -m fatou_backend.cli sexp --base e --values 0.5 --dps 100 --fatou-gp fork`;
  `FatouGP(n_workers=N)` parallelizes large batches.
- **Plain PARI/GP** (no Python, no state cache):
  ```
  \p 120
  \r src/fatou_backend/vendor/fatou_fork.gp
  sexpinit(2, 120, 4, 0);   \\ base, nlim ~ dps, nskip, looplim=0 (= converge fully)
  sexp(0.5)
  ```

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

- `fatou.gp` © Sheldon Levenstein ("sheldonison"), published on the
  Tetration Forum in the thread
  ["new fatou.gp program"](https://tetrationforum.org/showthread.php?tid=1017)
  (2015). It is vendored unmodified as `src/fatou_backend/vendor/fatou.gp`;
  the optimized fork `fatou_fork.gp` is a derivative work of it. The
  original carries no explicit license — see
  `src/fatou_backend/vendor/NOTICE.md` for the full attribution and terms
  status. This project exists thanks to that work.
- Base-change mathematics: the mixed-base tetration paper series
  (Janis Justus, 2026), included in `papers/` and permanently archived
  on Zenodo (DOIs in `papers/README.md`; © the author, all rights
  reserved).
- Everything else — wrapper, fork optimizations, references, atlas,
  tools — is © 2026 Janis Justus under the [MIT license](LICENSE)
  (the two vendor files and the papers are excluded, as stated there).
