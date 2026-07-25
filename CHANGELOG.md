# Changelog

## v0.1.0 — 2026-07-25

First tagged release. A Python wrapper and optimized fork of Sheldonison's
`fatou.gp` (Kneser tetration in PARI/GP), plus proven reference values and a
base atlas.

### Engine

- **30 gate-verified optimizations** in `src/fatou_backend/vendor/fatou_fork.gp`:
  persistent walk/Schröder caches, incremental Taylor/theta/extraction updates at
  reduced precision, FFT/Bluestein extraction on exact-size grids, radix-r mixed
  DFT, per-level degree truncation, precision laddering, contour-radius tuning.
- Every keep had to pass a **frozen accuracy gate** (`research/gate.py`, immutable
  together with `research/reference/`) — eleven checks over six real and complex
  bases plus three roundtrips. Thresholds are per-base and calibrated against
  the unmodified original, not uniformly "full digits"; the base-e groups sit at
  ~59 digits because that is what the engine truly delivers there.
- All DFT twiddle tables are built by doubling (~log₂(n) ulp) rather than by
  repeated multiplication (~n ulp). Measured against high-precision truth at
  precis 404, the doubling build is 3.3× better at n=192 and 34× at n=4608,
  with the advantage growing in n. `mixdft` caches its tables (151 hits / 10
  misses per run at dps 300), which is what makes the accuracy free — the
  change is clock-neutral to within 0.9% at dps 300 and 400.
- Two of the keeps remove *precision loss* rather than work, so the same `dps`
  now yields more true digits than before:

  | dps | true digits before | after |
  |---|---|---|
  | 300 | 295.0 | **302.2** |
  | 400 | 383.8 | **398.0** |
  | 520 | 496.9 | **508.9** |

- Measured speedup of the final round, paired back-to-back per row: **~2.3×** at
  dps 300/400/520, on top of the ~13× the fork already had over the unmodified
  original.
- Memory: the dps-1020 tier, which crashed the original engine, completes.

### Verification

- `research/reference/values.json` — reference values with *proven* error bounds
  via the error-vector method: `sexp_e(0.5)` to **972 digits**, plus 698/497
  tiers and `sexp_2(0.5)` to 497 (error vector **and** engine diversity).
- `research/tools/digits_vs_reference.py` measures true digits against those
  references and refuses to report agreement past the proven ceiling. Its
  `--pin` option pins the run to a fixed core set and verifies the mask took
  effect, refusing to report timings it could not actually pin.
- The DFT/twiddle layer is tested directly (`tests/test_transform_layer.py`):
  against a naive O(n²) DFT, against high-precision truth for the power
  tables, and by a source check that no `powers()` call survives at a twiddle
  site — the construction that carried ~n ulp instead of ~log₂(n).
- Tests: 72 passing, 4 skipped (`python -m pytest tests/`).

### Base atlas

- `fatou_backend.basechange` — tetration for any base from a single base-e
  anchor via a Fourier-mode table; ~1e-24 absolute with the default table,
  ~1e-37 with the k=40/dps-113 table, `slog(sexp(y))` round-trips to ~1e-52.
- Zero-dependency interactive demo:
  `research/tools/tetration_calculator.html`.

### Documentation

- `research/METHODS.md` — verification method, the full keep stack with
  mechanisms, measured speedups, and the negative results (what did *not* work
  and why), including the routes that were closed against the p^4.1 exponent.

### Known limitations

- **The depth-scaling exponent is ~4.1 and this release does not change it.**
  All optimizations here are constant-factor by construction; the local
  exponents confirm it (300→400: 3.31 before, 3.24 after). The structural reason
  is documented in `research/METHODS.md` §3: with Θ(p) Picard iterations each
  touching Θ(p²) digits, Θ(p³) is a floor for this discretization.
- The 520→1020 pair, where the 4.13 figure was originally measured, has **not**
  been re-measured for this release; an attempt was discarded as contaminated
  (see below).
- **Timings are sensitive to CPU core placement.** On a hybrid-core CPU
  (i7-12700H, 6 P + 8 E) the same engine and input measured 51 s pinned to
  P-cores vs 85 s left to the Windows scheduler — a 1.68× swing, digit-identical.
  Compare only paired, pinned, back-to-back runs. See the measurement-hygiene
  note in `research/METHODS.md` §2.
- The calibration rule is now **≥ 0.95·dps** true digits, measured at dps
  300/400/520. The older `dps − 24` rule is retired: it under-promises below dps
  500 and over-promises above it. The CLI's display floor deliberately stays at
  the conservative `dps − 24`.
- Bases `0 < b < 1` are not supported.

### Credits & license

`fatou.gp` © Sheldon Levenstein ("sheldonison"); vendored unmodified, and the
fork is a derivative work of it — see `src/fatou_backend/vendor/NOTICE.md`.
The papers in `papers/` are © the author, all rights reserved. Everything else
is MIT (`LICENSE`).
