# Methods — how the results in this repository were obtained and verified

This is the curated English record of the optimization and verification work
on the `fatou.gp` fork. It condenses a much longer internal lab journal;
everything stated here was measured on the frozen benchmark/gate of this
repository and is reproducible from the committed code.

## 1. Verification: how a digit claim is proven

No result in `research/reference/values.json` relies on "the engine says so".

**Error-vector method.** Kneser tetration via fatou.gp converges
iteratively; running the *same* computation twice at different iteration
depths and comparing the results bounds the true error of the worse run:
if two runs agree to `d` digits, the earlier run's true error is at most
~10^−d (the iteration is contractive; agreement cannot be accidental at
hundreds of digits). All "proven digits" figures use such pairs.

**Engine diversity.** For the 500-digit tiers the fork was additionally
checked against the *unmodified original* `fatou.gp` at full depth —
two different code paths agreeing to 495+ digits is the strongest
practical proof that the fork's 21 optimizations did not bend accuracy.

**Calibration law.** Across all measured tiers, true correct digits ≈
`dps − 24` (with 0–20 digits of quantization overshoot on top). Measured
convergence rates: base e ≈ 2.05 digits/iteration, base 2 ≈ 1.27.
This law is what the atlas certification prototype uses for its engine-error
term (e.g. dps 60 → engine ε ≈ 1e−33, *not* the naive 1e−55).

Every entry in `values.json` records its full provenance in the meta
block (engine, dps pair, measured error-vector delta, runtimes) —
German strings, but complete.

**The v2 reference correction** (referenced as `DECISION_reference_v2.md
Option A` in `gate.py` and the `values.json` meta; the original decision
record lives in the private lab notes): the initial v1 base-e reference
values beyond ~digit 64 were *correlated-error artifacts* — the wrapper
passed an iteration cap of 30 everywhere, base e converges at ~2.1
digits/iteration, so every v1 base-e number carried only ~64 true digits
while self-comparisons showed spurious full-precision agreement
(identical iteration trajectories on both sides). Decorrelated runs at
different iteration caps exposed this. With user approval ("Option A",
2026-07-10) the base-e entries were replaced by converged,
error-vector-verified values and the gate thresholds recalibrated, then
re-frozen. This episode is why the error-vector method (§ above) is the
only verification this project accepts.

**Proven reference ladder** (`research/reference/values.json`):

| value | proven digits | method |
|---|---|---|
| sexp_e(0.5), dps 1020/1033 pair | **972** | error vector |
| sexp_e(0.5), 700 tier | 698 | error vector |
| sexp_e(0.5), 500 tier | 497 | error vector + engine diversity (496) |
| sexp_2(0.5), 500 tier | 497 | error vector + engine diversity |

## 2. The fork: 21 gate-verified optimizations

Every change had to pass a frozen accuracy gate (all bases, full digits,
immutable reference values) before being kept. Measured end-to-end
(same machine, base e):

| target | true digits | original | fork | speedup |
|---|---|---|---|---|
| dps 220 | ~196 | 461 s | 73 s | 6.3× |
| dps 300 | ~276 | ~27 min | ~2.8 min | 9.5× |
| dps 520 | ~495 | ~5.4 h | ~24.8 min | 13× |
| dps 1020 | ~973 | memory crash | ~6.7 h | — |

Base 2: dps 300 went 54 min → 17.6 min.

The keeps, grouped by mechanism:

1. **Persistent caches across iterations.** The dominant cost is
   re-evaluating slowly-changing quantities from scratch each iteration:
   per-sample caches for the Schröder/inverse-Schröder walks
   (`isuperf`/`isuperf2` — the single largest win), the superfunction
   samples used by the theta map, and the Bluestein chirp vectors.
2. **Incremental updates at reduced precision.** Between consecutive
   iterations the Taylor coefficients, theta series and extracted
   coefficients change only in their low-order digits. The fork computes
   *differences* at a reduced working precision (`dps − re + margin`) and
   adds them to the cached full-precision values. Additive updates are
   cancellation-free, so no accuracy is lost — this is a scale-separated
   refinement scheme (the naive version, re-running at low dps and adding
   the correction, fails: the correction rounds away).
3. **Exact-size FFT extraction.** Coefficient extraction and the theta DFT
   run on exact-size grids via Bluestein's algorithm, with the underlying
   convolution done by power-of-two FFT; grid sizes are quantized (64/32
   steps, 512 above grid 3584) so caches stay valid across iterations.
   Unified for all bases, real and complex.
4. **Contour/sampling radius tuning.** The sampling radius for the
   family-of-circles contour was rescanned (`ctr·9/10` optimum); the stop
   criterion tightened accordingly (~6.2 stop-terms per digit instead of 9).
5. **Precision laddering.** Early iterations run at reduced precision with
   a guard that doubles back if the ladder would touch accuracy.

Total asymptotic cost is unchanged (~p^4.1 in the digit count p,
decomposing as iterations p^1.0 × grid² p^1.85 × arithmetic p^1.29);
the fork wins a large constant factor plus the memory fixes that make
dps 1020 possible at all.

One regression escaped the gate's coverage and was found later: the
"all bases" cache extensions (keeps 18–20) silently broke real bases
*below* e^(1/e) (attracting-fixed-point regime), which the gate never
tests. Fixed by routing those bases through the exact pre-extension
code path (`subeta` flag); verified against the unmodified original
engine (agreement ~1e-33 at b = 1.2) and guarded by a regression test.

## 3. Negative results (measured, not folklore)

Attempts to break the p^4.1 exponent, all benchmarked on this code base:

- **Aitken Δ² on coefficient sequences** — no acceleration (revert).
- **Anderson/Krylov acceleration of the fixed-point iteration** — the
  iteration *is* affine on a fixed grid (verified to machine precision),
  but its spectrum is uniform (~1.67–1.74 digits/step across all modes),
  so Krylov needs k ~ p dimensions: 1.9× fewer steps at best, no
  exponent change.
- **Fast multipoint evaluation (FLINT)** of the orbit walks — numerically
  unstable on the disk-spread walk endpoints (product-tree factors ~0.55
  per level collapse to underflow at depth 2048); parked pending
  shell-batching.
- **Nested/tripled grids** — offsets nest exactly, but the 3× grid
  overshoot costs more than the reuse saves.
- **Tail truncation of the difference polynomials** — the differences are
  full-band; only 2.2% savings available.

These closures are load-bearing: they mean the remaining routes to a lower
exponent are algorithmic restructuring (batched multipoint) or the
theory-level approaches of the paper series, not tuning.

## 4. The base atlas

`fatou_backend.basechange` evaluates sexp/slog for *any* base from a single
base-e engine, using the base-change phase function Φ of the paper series
(see `papers/`):

- **Peel ladder** for Φ_{e,b}(θ): tower up explicitly to >1e4, enter the
  target-base world by the exact two-level identity
  W = log_b(log_b T) = α̂ + c·y (no rounding of the huge tower),
  peel down numerically, one slog call.
- **Mode table** (`research/reference/phi_modes.json`): µ plus a k ≤ 20
  Fourier-mode head per base at dps 60, computed on demand (~0.2 s) and
  cached (~3 KB/base). Reconstruction error vs. direct Φ: <1e−20.
- **Anchor paths**: `sexp_anchor` lifts the height by n + Φ and descends
  in ln-coordinates (truncation is sub-precision above the threshold);
  `slog_anchor` climbs b-towers, enters the e-world exactly, and solves
  x + Φ(x) = target by Newton (Φ′ small, 3 steps suffice).
- **Validation** against the proven references: ~1e−24 absolute with the
  default table (k=20, dps 60), ~1e−37 with k=40/dps 113; round-trips
  slog(sexp(y)) agree to ~1e−52. The µ-hub
  (`research/reference/mu_hub.json`, 26 bases in [1.5, 100]) matches
  published anchor values to 31–32 digits.
- **Certified enclosures** (prototype, `research/tools/m54_cert_proto.py`):
  ball arithmetic over the grid with three explicit error terms — engine
  error from the calibration law (§1), grid aliasing (measured exactly:
  A_16 ≈ 1.3e−20), and interval width.

## 5. Reproducing

- **Any reference value**: `python research/tools/verify_reference.py
  --key "sexp|e|0.5|500" --dps 520 --pair 533` re-runs the engine and
  reports both agreement with the stored value and a fresh run-vs-run
  error vector. Low tiers take seconds-to-minutes; the deep tiers cost
  what the table in §2 says.
- `bench/` — frozen workload and metrics; the headline benchmark
  (`bench/results/2026-07-12-m8-keeps21.json`) references the vendored
  fork by relative path and reruns as-is.
- `research/gate.py` + `research/reference/` — the frozen accuracy gate
  used for every keep decision (gate and references are immutable).
- `tests/` — includes atlas validation against the proven references
  (`FATOU_BACKEND_RUN_SLOW=1` enables the slow gate block).
- The atlas artifacts (`phi_modes.json`, `mu_hub.json`,
  `research/tools/calc_data.json`) regenerate from scratch via
  `basechange.phi_modes_cached`, `research/tools/m53_mu_hub.py` and
  `research/tools/calc_export.py`.
