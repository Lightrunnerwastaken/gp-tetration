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

**Calibration law.** True correct digits ≈ `precis − (precis − 4.7)/21`,
where `precis` is PARI's actual working precision (word-quantized: dps 300
gives 308, dps 1020 gives 1021). Measured against the proven references:

| dps | precis | measured | this law | the older `dps − 24` |
|---|---|---|---|---|
| 300 | 308 | 295.0 | 293.6 | 276.0 |
| 400 | 404 | 383.8 | 385.0 | 376.0 |
| 520 | 520 | 496.9 | 495.5 | 496.0 |
| 1020 | 1021 | 973.0 | 972.6 | **996.0** |

The loss is **proportional**, not constant: `precis/21`, from a cancellation in
the Schröder walk. `isuperf` iterates to within `isuperfr = 10^(−precis/21)` of
the fixed point and then evaluates `subst(fsl, x, y − L)` — `y` and `L` are both
O(1) while their difference is at scale 1e−19 at dps 400, so an absolute
1e−precis error re-enters as a *relative* 1e−(precis−19) one. The 21 is the
file's `\ps 21` series precision.

The earlier `dps − 24` form is conservative below ~dps 500 — which is why the
atlas certification, which uses it at dps 60, stays safe (36 claimed against
~74 real) — but it **over-claims by 23 digits at dps 1020** and must not be
used to size a deep run.

Measured convergence rates: base e ≈ 2.05 digits/iteration, base 2 ≈ 1.27.

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

## 2. The fork: 27 gate-verified optimizations

Every change had to pass a frozen accuracy gate (all bases, full digits,
immutable reference values) before being kept. Measured end-to-end
(same machine, base e):

| target | true digits | original | fork | speedup |
|---|---|---|---|---|
| dps 220 | ~196 | 461 s | 73 s | 6.3× |
| dps 300 | ~276 | ~27 min | ~2.8 min | 9.5× |
| dps 520 | ~495 | ~5.4 h | ~24.8 min | 13× |
| dps 1020 | ~973 | memory crash | ~6.7 h | — |

Base 2: dps 300 went 54 min -> 17.6 min.

Absolute times in that table are only comparable *within* one measurement
session: the same engine that shows ~24.8 min at dps 520 above measures
9.2 min on a later machine. Only the ratios travel.

Six further keeps (2026-07-25) add a factor that grows with depth, measured
against the proven references on one machine in one session:

| dps | before | after | speedup | true digits before/after |
|---|---|---|---|---|
| 300 | 83.95 s | 40.28 s | 2.084x | 295.0 / 295.2 |
| 400 | 216.38 s | 97.73 s | 2.214x | 383.8 / 387.9 |
| 520 | 555.05 s | 291.00 s | 1.907x | 496.9 / 496.7 |

The digit column moves because the last of these keeps removed a cancellation
rather than an operation (see section 1): at dps 400 the same run now carries
4 more true digits than before, so the speedup at a fixed *digit* target is
larger than the table's fixed-dps ratio.

The 0.4-digit gap at dps 520 (verified against the 972-digit reference, not the
497-digit one) is where the two trajectories stop, not a loss: both sit above
the calibrated floor of dps - 24 = 496, inside the 0-20 digit overshoot band
of section 1, and the digit counts are identical at the other three tiers.

They are constants, not an exponent change, and the locally measured exponents
say so: 3.29 -> 3.04 for the 300->400 pair and 3.61 -> 3.91 for 400->520. The
movement in both directions is grid-quantization noise; nothing shifted the
exponent, which is what section 3 argues must be the case.

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
6. **Decision-only evaluations are not evaluations.** In the sampling loop the
   theta branch used a full O(N) Horner solely to settle three *discrete*
   decisions (which theta map, and two period-representative comparisons).
   The measured margins are O(1) -- |Im| >= 0.795, relative period-tie margin
   >= 0.9998, and the period loops shift zero times at any precision -- so one
   evaluation per grid stretch reproduces every decision bitwise.
7. **The theta extraction had been left behind.** thtaylor still ran the
   quadratic rotation DFT that the main extraction lost early on; a section
   profile (rebuilt against the current engine, not an old one) put it at
   29.5% of the run at dps 200. Same transform, FFT/Bluestein.
8. **Radius re-tuning after the balance moved.** The sampling radius optimum
   had been scanned once and confirmed once; ten keeps later the N^2 sampling
   term carries a larger share and the optimum had moved. Re-scanning is
   cheap (the knob already exists) and gave a factor that grows with depth.
9. **Evaluate only the degree each point needs.** The coefficients decay like
   circr^-k, so a sample at radius rho needs ~dig/log10(circr/rho) terms --
   but every sample was getting the full degree, and the samples span
   |w| = 0.07 .. 0.96. Bucketing them by radius (the bucket is a property of
   the cached walk endpoint, so it is derived once per grid stretch) and
   keeping one truncated copy of the polynomial per bucket cuts the Horner
   work 2.4x on a measured state. The truncation index comes from the actual
   coefficient magnitudes, not from a decay model.
10. **A quantization guard that never fired where it mattered.** The theta
   grid is quantized so its caches can engage -- but only above 64 samples,
   and the run spends a contiguous block of iterations below that, where the
   grid moved every pass, the caches were reallocated and everything was
   recomputed at full precision. Measured: 60 of 232 iterations at dps 300.
11. **The grids are never arbitrary, so Bluestein is overkill.** Exact-length
   DFTs went through Bluestein, which needs FFTs four to eight times the data
   length. But the grid quantum is a power of two, so every length factors as
   r*2^k with an odd r <= 9 -- a radix-r decimation gives r power-of-two FFTs
   of length N/r plus N*r twiddles instead (~82k complex multiplications
   against ~688k at the largest size). It is an exact reorganization of the
   same sum, verified to agree with Bluestein to 1e-70 before being timed.

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
  unstable on the walk endpoints. First read as an underflow of the
  product-tree factors (~0.55 per level) and parked pending shell-batching;
  re-measured in 2026-07-25 with shells and honest guard budgets, which
  closed it for a different and stronger reason — see the round below.
- **Nested/tripled grids** — offsets nest exactly, but the 3× grid
  overshoot costs more than the reuse saves.
- **Tail truncation of the difference polynomials** — the differences are
  full-band; only 2.2% savings available.

Later round (2026-07-25), aimed squarely at the exponent:

- **Product-tree multipoint, honestly budgeted** -- the classical fast
  evaluation, tested on a real frozen engine state rather than synthetic
  points. It needs about **6 guard bits per point**. The geometry sets the
  amplification -- the evaluation points are the images of ONE circle under
  three analytic maps (identity, exp, log), so they sit on smooth arcs, the
  worst case for subproduct polynomials -- and the guard budget then buys it
  back exactly: on a frozen dps-300 state with 1280 points, 4.0*N guard bits
  still returns garbage while 6.0*N gives bitwise agreement with Horner.
  Because that demand is linear in N, the precision growth cancels the
  algebraic gain, and with an honest guard the tree measures several times
  *slower* than Horner across the working range.

  Worth recording, because it misdirected the follow-up work: this was first
  reported here as "a guard sweep from 0 to 20480 bits changes nothing, so the
  cause is geometry rather than budget". That was an artifact of the measuring
  harness, which built the ball inputs once at a fixed precision and then
  raised only the working precision -- the frozen input radii pinned the result
  no matter how much guard was added. Rebuilding the inputs at each precision
  gives the law above. The route stays closed, on timing rather than on
  stability.
- **Constant-precision increment ladder** -- the incremental path evaluates the
  diff at `working precision - residual + 40` digits, which is linear in the
  residual because the working-precision ladder is `2*re + 60`. A factor of 1
  would make it constant and, in the cost model, move the exponent from ~3.9
  to ~3.05. It measures 2.19x faster at dps 300 and costs 82 true digits.
  Achievable digits saturate at ~115*m + g + 38 with a *fixed* offset, so the
  margin needed grows with the target precision and the gain vanishes. The
  doubling is load-bearing: the working precision has to reach full precision
  early enough (halfway, at re = (p-60)/2) for the state to resolve at all.
- **Universal asymptotic form for the theta harmonics** -- the one route that
  attacks the iteration count rather than the cost per iteration. A best-fit
  three-parameter form reproduces log10|theta_m| to ~1e-3 decimals
  out-of-sample and the residual *grows* with m; more parameters are measurably
  worse out-of-sample. A later parameter-free reading of the same data is
  sharper and reaches the same verdict: the m-exponent measures 1.000004
  (0.999999..1.000054 over m = 40..115) and the amplitude equals the engine's
  own constant rlnlm2 = 1/L2 to 9e-7 -- so the law is simply the Darboux term
  of the logarithm the construction already writes down,
  `theta_m = -(1/L2) * w*^-m / m` with `|w*| = 257.81`. That is 6.76 predicted
  digits per harmonic, i.e. seeding it saves ~2.8 harmonics out of ~0.78p:
  a constant fraction of a percent of the iterations, constant in p.
- **Lowering the direct-vs-theta branch radius** -- moves samples onto the
  branch that (after keep 6 above) needs no polynomial evaluation at all. It
  gives 1.14x and 22% fewer iterations at dps 200 and *reverses* to 0.94x at
  dps 300, because the theta series converges more slowly closer to the
  centre. Pushing further breaks accuracy outright (92 instead of 193 digits).

- **Series composition for the exp arc** -- the evaluation points turn out to be
  the images of ONE circle under three analytic maps (identity, exp, log), so
  the exp arc's values are in principle one series composition plus one FFT,
  in `O(M(N) log N)` and with no C dependency (PARI's own `fft` suffices).
  Measured before building it: on the sampling circle the composed function
  already reaches **10^1141** while the values actually used on the arc are
  O(1), and it needs 4.22*N terms rather than the estimated 2.9*N.

**The mechanism behind all of these.** The evaluation points lie on *arcs*, not
on the full circle. Every global representation -- a product tree over all
points, a series composition over the whole circle, a low-rank/FMM compression
-- has to represent the function where it is astronomically large and then
cancel back down to O(1) on the arc. The price is always a guard or rank
requirement proportional to N, and that consumes the algebraic gain exactly.
A viable candidate would have to work *only on the arc* and still be
quasi-linear; the only one we know of is Moroz (FOCS 2021), and whether it
tolerates points on arcs is open.

These closures are load-bearing: they mean the remaining routes to a lower
exponent are algorithmic restructuring or the theory-level approaches of the
paper series, not tuning.

**Where the exponent actually stands.** The iterate carries N(p)*p ~ 6.2 p^2
digits, and the iteration count is Theta(p) -- measured, and closed from four
independent directions (extrapolation, Krylov, spectral, asymptotic form). Any
scheme that touches its state once per iteration therefore costs Theta(p^3).
In the cost model an N^2 evaluation gives exponent 3.8-3.9 (measured 4.13 on
the 520->1020 pair) and a hypothetical stable quasi-linear evaluator would give
2.8-2.9. So ~2.9-3.3 is the reachable band for this discretization, and getting
there needs an evaluation primitive that survives points on analytic arcs --
the classical product tree demonstrably does not.

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
