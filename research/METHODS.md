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
practical proof that the fork's optimizations did not bend accuracy. Note the
scope: that comparison was run against the 10-keep engine of the time, so it
does not by itself cover the later keeps. What covers those is the frozen gate
after each one, plus `tests/test_fork_engine.py`, which re-runs the diversity
comparison (bases e, 2, 1+I) against the current fork on every test run.
Base 1.2 was in that list until 2026-07-30 and has been removed: below eta
the two engines share a construction, so their agreement measured the
construction rather than the answer. It is now checked against an
independent implementation instead (§2).

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

**The law above describes the engine up to exp-056.** Two later keeps changed
it: exp-057 removed the cancellation itself (offset iteration in `isuperf`),
and exp-060 replaced the file's fixed `\ps 21` by
`myps = max(seriesprecision, min(64, round(precis/8)))`, so the denominator is
no longer a constant. The closed form has **not** been re-derived for the
shipped engine; use the measured table:

| dps | true digits, exp-056 engine | true digits, shipped engine |
|---|---|---|
| 300 | 295.0 | **302.2** |
| 400 | 383.8 | **398.0** |
| 520 | 496.9 | **508.9** |
| 1020 | 973.0 | ≥972 (reference exhausted) |

Conservative rule for sizing a run: **true digits ≥ 0.95 · dps**, verified at
dps 300/400/520. At dps 1020 the run saturates the 972-digit proven reference,
so the calibration there is not measurable with the current ladder — the entry
is a floor, not a measurement. `research/tools/digits_vs_reference.py` enforces
that ceiling and reports `>=972 (reference exhausted)` rather than the raw
agreement, which is the number an unbounded comparison would print.

Measured convergence rate, base e: **1.30 digits/iteration** (dps 100: 83
iterations to 110.96 digits; dps 200: 158 to 204.92; marginal 75 iterations for
93.96 digits). The earlier figure of 2.05 published here described the
pre-exp-023 engine and was 58% too high — exp-023 and exp-045 traded convergence
rate for a smaller grid, which is the whole point of those keeps, and the
published rate was never updated. In the cost law this makes the iteration count
**I(p) ≈ 0.80·p**, not 0.59·p. The *exponent* is unaffected — I(p) is Θ(p) either
way — but any absolute cost estimate built on 0.59 or on 2.05 is wrong.

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
| sexp_e(0.5), dps 2000 run | **1973** | error vector, 2000/2033 pair (2026-07-30) |

The last row is a different kind of evidence and is worth separating from the
error-vector rows above it. A single dps-2000 run (26.69 h CPU, polynomial
degree 10752, `research/candidates/candidate_e_2000.txt`) reproduces every
entry of the ladder *at or above its certified depth* — 509 / 717 / 992
digits against the 500 / 700 / 1000 tiers whose proven depths are 497 / 698 /
972. It also settles the previously unproven `candidate_e_1020`: the two
agree to **973** digits, and that candidate's own contour residual claimed
973.3. Two computations at different precision, on grids of 7168 and 10752
terms, landing on the same digit count the weaker one predicted for itself is
the strongest cross-check the project has produced so far.

That comparison could not establish the depth of the dps-2000 value itself —
an agreement is capped by the weaker partner, and nothing deeper than the
1000 tier existed. The dps-2033 partner has since been computed (20.60 h CPU,
degree 11264, contour residual 2010.8): the pair agrees to **1973** digits,
98.7 % of dps. Two controls keep that from being self-confirmation. Both runs
agree with the proven ladder *identically* — 509 / 717 / 992 against the
500 / 700 / 1000 tiers — so they do not share an error at those depths; and
the 2033 run's own contour residual (2010.8, i.e. 98.9 % of its dps) lands
where the pair predicts. The earlier extrapolation from re/dps = 0.954
suggested ~1908 and was conservative.

One caution the pair also produced: the two runs took 26.69 h and 20.60 h
for essentially the same computation — the *higher* precision was 23 % faster.
The logs differ in stack behaviour (the slower run regrew small stacks
repeatedly; the faster one climbed once to 1 GB), but that is a partial
explanation at best. Treat absolute runtimes in this document as accurate to
no better than ±25 % unless a drift bracket is quoted with them.

## 2. The fork: 30 gate-verified optimizations

Every change had to pass a frozen accuracy gate (eleven checks over six real
and complex bases plus three roundtrips, against immutable reference values)
before being kept. The per-base thresholds were calibrated once against the
unmodified original, so they are *not* uniformly "full digits": the base-e
groups require 58.9 and 59.2 digits against 80- and 200-digit references,
because ~64 digits is what the engine actually delivers for base e at those
settings after the v2 reference correction. Base 2 at dps 200 requires 189.

Measured end-to-end (same machine, base e):

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

**Measurement hygiene (measured, not assumed).** On a hybrid-core CPU the
comparison can break even *within* a session. On an i7-12700H (6 P-cores +
8 E-cores) the same engine and input measured

| run | dps 300 |
|---|---|
| pinned to P-cores (`ProcessorAffinity = 0xFFF`) | 51 s |
| left to the Windows scheduler | 85 s |

— a **1.68× swing from core placement alone**, digit-identical, and it appears
when foreground applications start and push a long-running background job onto
the E-cores. A further ~1.4× drift accumulated over a day of sustained load
(thermal), so a morning number and an evening number of the same engine differ
by ~2.3×. Both effects cancel in a *paired* ratio only if both runs are pinned
and back-to-back; a single unpinned run at one precision compared against an
earlier run at another precision measures the scheduler, not the algorithm.
Every speedup ratio quoted here is from paired runs; the depth-scaling exponent
is quoted only from pinned pairs.

Nine further keeps (2026-07-25) add a roughly constant ~2.3×, measured against
the proven references, each row a paired back-to-back run of both engines:

| dps | before (b0ecd1f) | after | speedup | true digits before/after |
|---|---|---|---|---|
| 300 | 84.64 s | 36.25 s | 2.34× | 295.0 / **302.2** |
| 400 | 219.31 s | 92.62 s | 2.37× | 383.8 / **398.0** |
| 520 | 554.00 s | 242.84 s | 2.28× | 496.9 / **508.9** |

Re-verified later the same day on a slower machine state, both engines pinned
to P-cores — the absolute times move by ~1.4×, the ratio does not:

| dps | before | after | speedup | true digits |
|---|---|---|---|---|
| 300 | 116.28 s | 51.47 s | 2.26× | 295.0 / 302.2 |
| 400 | 301.67 s | 130.64 s | 2.31× | 383.8 / 398.0 |

The digit column moves because two of these keeps removed precision losses
rather than operations (exp-057, exp-060; see section 1), so the gain at a
fixed *digit* target is larger again than the fixed-dps ratio: the dps-520 run
that used to deliver 497 digits now delivers 509.

The 0.4-digit gap at dps 520 (verified against the 972-digit reference, not the
497-digit one) is where the two trajectories stop, not a loss: the digit counts
are identical at the other three tiers, and both sit above the measured floor
for that tier. (This paragraph previously justified the gap against
`dps - 24 = 496`; that rule is retired -- see section 1.)

They are constants, not an exponent change, and the locally measured exponents
say so. From the pinned pair above, 300->400: **3.31 before, 3.24 after** —
and the unpinned morning pair independently gives 3.31 -> 3.26. The 400->520
pair moves the other way (3.53 -> 3.67). The movement in both directions is
grid-quantization noise; nothing shifted the exponent, which is what section 3
argues must be the case.

**The 520->1020 pair, re-measured -- twice, and the second one supersedes the
first.** Attempt one (unpinned, loaded machine) was discarded outright.
Attempt two was pinned and bracketed -- dps 520 immediately before and after the
1020 leg -- and gave 4.35 in an interval [4.23, 4.48], the width coming from an
18.7% slowdown the bracket caught across that run.

Attempt three was not aimed at the exponent at all: it is the dps 1200/1265
error-vector pair (section 1), which happens to span the same depth range on a
machine that stayed still. Bracket drift 0.991 -- i.e. 0.9% *faster* at the end.
Three intervals from that one session:

| interval | exponent |
|---|---|
| 520 -> 1200 | 4.159 |
| 520 -> 1265 | 4.151 |
| 1200 -> 1265 | 4.025 |

**These supersede the 4.35.** The contradiction is decisive rather than a
matter of taste: the exponent grows with depth, so 520->1200 must exceed
520->1020. Measured, it is lower (4.16 against 4.35). What differed is the
hardware state, not the algorithm -- the same dps-520 workload took 348-413 s
during attempt two and 243 s during attempt three, so that machine was
1.4-1.7x slower throughout. A bracket catches drift *within* a run; it cannot
tell you the whole session is running degraded, and attempt two's was.

So: **the exponent is ~4.0-4.16, unchanged from the historical 4.13**, which is
what every keep being constant-factor by construction predicts. The local
exponent still grows with depth: 3.24 (300->400), 3.39 (300->520), ~4.16
(520->1200).

The lesson is about method rather than about tetration. Two hours of sustained
load moved this machine by 18.7%, and a whole session can sit 1.4-1.7x off its
own best state. Absolute timings compare only between runs minutes apart, and
an exponent is only worth quoting when a bracket shows the machine held still --
0.991 here, against 1.187 for the number this replaces.

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
11. **A cancellation, not an operation.** The Schroeder walk ran to within
   10^(-precis/21) of the fixed point and only then formed `y - L`, turning an
   absolute error into a relative one and costing precis/21 digits -- which was
   the whole of the published calibration law. Iterating the offset `u = y - L`
   with expm1/log1p removes it: same algorithm, different coordinate, and the
   run comes out both faster and four digits more accurate at dps 400.
12. **A series length frozen at a constant while the precision grew.** The
   Schroeder walk runs until |u| <= isuperfr, and log10(1/isuperfr) is exactly
   precis/seriesprecision -- so a file-level `\ps 21` pins the walk at
   precis/21 decades and, with it, a precision loss of the same order. Scaling
   the series with the working precision buys those digits back and the
   shorter walk pays for the longer series build: across ps = 21..64 the clock
   is flat while the true digits climb from 387.9 to 400.2.
13. **A geometric sequence built with N exponentials.** Both sampling grids are
   `c * mu^s` (x1 is affine in s) but were built point by point with a
   full-precision complex exp, and thtaylor then undid each one with a log.
   A doubling build keeps the ulp growth at O(log N) instead of O(N).
14. **The grids are never arbitrary, so Bluestein is overkill.** Exact-length
   DFTs went through Bluestein, which needs FFTs four to eight times the data
   length. But the grid quantum is a power of two, so every length factors as
   r*2^k with an odd r <= 9 -- a radix-r decimation gives r power-of-two FFTs
   of length N/r plus N*r twiddles instead (~82k complex multiplications
   against ~688k at the largest size). It is an exact reorganization of the
   same sum, verified to agree with Bluestein to 1e-70 before being timed.
   **That verification was insufficient** -- see below.

**A correctness fix, not a speedup (exp-061).** Every twiddle/root table in
the fork was built with PARI's `powers()`, i.e. by repeated multiplication,
which accumulates ~n ulp. Seven call sites did this, four of them on the main
power-of-two path, including one (`pw` in `thtaylor`) that is not an FFT input
at all but multiplies the output coefficients one by one. A doubling build
(`geopow`) produces the same vector with ~log2(n) ulp; measured against
high-precision truth at precis 404 it is 3.3× more accurate at n=192 and 34×
at n=4608, the advantage growing in n as n vs log2(n) predicts.

Why item 14's check could not see it: it compared the new radix-r DFT against
Bluestein, and *both* built their twiddles with `powers()`. The two agreed to
1e-70 on a common-mode error, and that agreement was read as proof. A
verification whose reference shares the suspect component verifies nothing.

The fix is free rather than costly because `mixdft` now caches its tables
(they depend only on grid size and precision, but were rebuilt every call):
151 cache hits against 10 misses per run at dps 300. Measured paired and
pinned, best of two: dps 300 46.97 s → 46.47 s, dps 400 120.45 s → 120.25 s,
true digits unchanged at 302.2 and 398.0.

The transform layer now has direct tests (`tests/test_transform_layer.py`):
against a naive O(n²) DFT, against high-precision truth, and — the guard that
actually catches this defect class — a source check that no live `powers()`
call remains. The numeric tests cannot catch it: at their dps-60 working
precision a 288-ulp error is ~3e-65 against a 1e-45 tolerance.

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

**That verification was weaker than it looked** (measured 2026-07-29). It
compared two implementations of the *same* construction, so it could only
confirm that they share an error, not that either is right. Against an
independent reference — regular iteration at the real attracting fixed point,
Koenigs/Schröder, no contour and no theta function
(`research/tools/regular_subeta.py`, self-tests S(1)=b and S(2)=b^b holding
to 148 digits) — both engines deliver about **15 correct digits at b = 1.2 when
60 are requested**, and nothing usable from dps ~100 upward. The tell is an
imaginary part of 7.26e−14 on a quantity that is provably real for
1 < b < e^(1/e). At dps 100 the two engines no longer even agree with each
other (degree 932 vs 1136), so the `subeta` isolation is incomplete above
dps 60 as well.

This is not a fork defect — the unmodified `fatou.gp` behaves the same way,
and sub-eta support was never its purpose. It is a *silent* one, which puts
it in the same category as the second defect below: a plausible number
returned with no error and no warning.

Two consequences worth stating plainly. First, the regime is cheap by the
right method and expensive by this one: regular iteration reaches 2000
verified digits for b = 1.2 in **100 seconds** (exponent ~2.6–3.0), against
26.69 h for base e through the Kneser machinery. Second, sub-eta bases cost
*more* here despite converging at the same rate (~1.25 vs ~1.30 digits per
iteration) because their sampling grid is twice as large at equal iteration
count, and the grid enters quadratically.

### Both defects fixed (exp-071, 2026-07-29)

**Sub-eta now uses regular iteration** (`regfix`/`regsigma`/`regsigmainv`/
`reginit`/`regsexp`/`regslog`; `sexp` and `slog` route there when `subeta`).
Measured against the independent reference: **403 of 403 digits** at a
400-digit target, with S(0)=1, S(1)=b and slog∘sexp all holding past the
target and a zero imaginary part. Cost at b = 1.2 fell from 20.5 s to 15 ms
at dps 60, and dps 200 went from "no result in 600 s" to 47 ms.

The entry point is `sexpinit(b)`, not `loop(kc)`: recovering b from
kc = 1 + log(log b) rounds the base to the ambient precision, and regular
iteration needs it to roughly *twice* the target (σ⁻¹ amplifies every
residual by λ⁻ⁿ). For the same reason the base must be exact — `sexpinit(1.2)`
now refuses with that explanation, and the Python backend converts decimal
sub-eta bases to rationals, scoped to that interval so no gate base changes.

**The second defect was `safefs`.** It is a search aid, not an evaluator:
above 5E8 it replaces the imaginary part with `random()`, and above
|Re z| > 10000 it returns the hardcoded sentinel `1E400*exp(I*imag(z))`.
Inside `abel()` that is correct and load-bearing — the walk only needs to
know "far away". But `invabel()` used it to *build* its return value, so the
sentinel became the answer: sexp_2(5) came back as 1E400/ln 2 =
1.4426950408889634E400 instead of 2^65536 ≈ 1e19728, and sexp_2(5.5) came
back negative. The reconstruction loop now uses the exact `fs()`; sexp_2(5)
is exact, and arguments genuinely beyond the representable range raise
`e_OVERFLOW` instead of returning a plausible number.

Both fixes are disjoint from the optimized path by construction, and the
frozen gate confirms it the strong way: 11/11 PASS with values *bit-identical*
to the pre-fix run (2.56296e−87 / 1.05895e−88 / 9.36335e−96), not merely
inside tolerance.

### What the sub-eta regime costs, and what it yields

Error-vector pairs (p, p+33) at b = 6/5, the same protocol as the base-e
ladder:

| pair | agreeing digits | share of dps |
|---|---|---|
| 250 / 283 | 250 | 100.0 % |
| 500 / 533 | 499 | 99.8 % |
| 1000 / 1033 | 999 | 99.9 % |
| 2000 / 2033 | **2002** | 100.1 % |

Against the independent mpmath implementation: 499 digits at dps 500 and
**2002** at dps 2000 — the two verifications agree *exactly*. That is worth
more than confirming the value: it calibrates the error-vector method itself
against a foreign computation, which the base-e ladder can never do for lack
of an independent implementation at depth. Where Kneser yields ~95 % of the
requested digits (973.3 of 1020), regular iteration loses essentially none.

### A third silent defect, and the pattern behind it (exp-074, 2026-07-30)

The band *just above* eta was broken too, and by a fork keep rather than
upstream. At b = 1.4494 (Re(Period) = 47.23) the engine returned **1.837**
contour digits where 150 were requested, printing 38 of them without complaint.
The unmodified original delivers 60.67 there.

Located by sweeping b = 1.4494 against all 34 fork versions at a uniform 240 s
budget: the collapse is a single commit, `af28712`, and an A/B switch inside it
isolated **exp-048** (per-point degree truncation) from exp-051. Turning the
truncation off restored 140.3 digits.

The cause is one anchor. exp-048 keeps term k while

    exponent(c_k) + k·log2(rho) ≥ mx − dig·log2(10) − 16,   mx = max_k exponent(c_k)

so the tolerance hangs on the largest *coefficient*, where it should hang on
the largest *contribution* |c_k|·rho^k. The two coincide exactly when the
coefficients decay — then max(lg2) is at k = 1 and (k−1)·log2(rho) = 0 — and
the function's own comment states that assumption: *"ct's coefficients decay
like circr^-k"*, which holds for **circr > 1**. Base e has circr = 1.3372;
this base has ≈ 0.133, so its coefficients *grow*, mx sat at index 160 of 161
while the largest contribution sat at index 2, and the threshold came out
284 bits ≈ 85 decimal digits too high. The kept degree fell from 161 to 45,
`sfunc` returned nonsense, and the residual stalled at once.

The fix anchors the threshold per level to that level's largest contribution.
Since bv ≤ mx always, it can only *lower* the threshold and keep *more* terms:
it cannot reduce accuracy anywhere. Measured at b = 1.4494, dps 150: 1.837 →
**147.410** digits, and the run now *converges* (looplim ≈ 143.5) instead of
stalling, with a value bit-identical to the unmodified original across all of
the original's valid digits. The gate is unchanged to the last digit, verified
both on a test copy via `gate.py --fork` and on the fork itself.

**It is not free, and the reason is not what it first appeared to be.** An
earlier draft of this section claimed the fix was a no-op wherever the old
assumption held, on the argument that decaying coefficients put max(lg2) at
k = 1, where (k−1)·log2(rho) = 0 and bv = mx. Measurement refuted that: for
base e the maximum sits at index 32, then 64, of 65 — **its coefficients grow
here too**. The old threshold was therefore already too aggressive for base e
(bv = −24.2 against mx = −7, some five decimal digits), merely without
consequence: the gate values are bit-identical because the extra terms do not
reach gate precision.

The cost is those extra terms, and it is **+9.7 % on base e** — median of five
alternating runs per variant at dps 150, ranges disjoint (5500–5890 ms against
6016–6516 ms) despite 7–8 % within-variant spread. A first attempt to buy the
correctness back by vectorising the search (`vecmax(lg2 + lr2·icsq)` instead of
an interpreted loop) changed nothing measurable: 9.7 % against 9.6 %. The
expense is the retained terms in the N² path, not the search.

Kept anyway, on the user's explicit decision (2026-07-30): correctness in a
regime the old code destroyed outweighs 9.7 % on one base. The rejected
alternative was a cutoff — use the cheap threshold while mx − bv stays small —
which would have restored the speed at the price of a magic constant encoding
exactly the kind of unstated assumption that produced this defect and the
sub-eta one before it.

An aside worth keeping: exp-048 had accidentally *masked* a non-termination
introduced by exp-018, which lifts the iteration cap until the digit goal is
met — for a base that never meets it, forever. Breaking the computation made
the loop exit. A stall that looks like a fix.

**The pattern, now twice confirmed.** exp-032 broke the sub-eta regime;
exp-048 broke the near-boundary band. Both were framed as "all bases", both
carried an unstated geometric assumption (the fixed point's nature; circr > 1),
and neither was visible to the gate, whose six bases all sit comfortably far
from eta. The question to ask of any future keep is therefore not only *does
the gate pass* but **what does this assume about the geometry, and which base
would violate it?** Both regimes now have their own regression test, since the
gate by construction cannot cover them.

Cost scales with an exponent near **2.5** rather than 4.1 — measured over
dps 250…4000, but with a **broken drift bracket** (the same dps-250 run took
46 ms at the start and 78 ms at the end, and dps 2033 came out *faster* than
dps 2000), so the number is an order of magnitude, not a measurement, and the
clean re-run is still owed. Structurally the lower exponent is what the method
predicts: regn ~ p steps, each one exp/log at working precision ~2p, so
p·M(2p) — with no sampling grid, no theta series, and hence no N² term.

## 3. Negative results (measured, not folklore)

Attempts to break the p^4.1 exponent, all benchmarked on this code base:

- **Aitken Δ² on coefficient sequences** — no acceleration (revert).
- **Anderson/Krylov acceleration of the fixed-point iteration** — the
  iteration *is* affine on a fixed grid (verified to machine precision),
  so Anderson is exactly GMRES on (I−A)x = b. Conclusion: **~2× fewer
  steps at best, no exponent change.** Re-measured end to end in
  2026-07-28 and confirmed, with a sharper picture than the first pass —
  see below.

  The first pass reported the spectrum as *uniform* (~1.67–1.74
  digits/step across all modes). Arnoldi on the scaled operator says
  otherwise: it is strongly **clustered** — |λ| = 0.051, 0.033 (×2),
  ~0.010 (×2), 0.0031 (×2), 0.00066, then a flat cluster at ~2.1e−4 —
  and scale-invariant (quadrupling the dimension adds no large
  eigenvalues). The check that settles it: the engine's own convergence
  rate falls *out* of the spectrum rather than being put in, four times
  independently (1.30 measured; 1.3113 and 1.2936 from the 128- and
  512-dimensional spectra; 1.2958 from a direct prototype run).

  Because the spectrum is clustered, GMRES converges superlinearly:
  digits(k) = 0.156 k² + 2.74 k, with the quadratic coefficient invariant
  under a 4× change of dimension, i.e. k ~ √(digits). In the running
  engine, with the grid ladder made geometric and GMRES run to each
  level's capacity, that reproduces as **k ≈ 2.3 √budget over six levels**
  spanning a factor 37 in budget.

  That is a genuine √p law for the *iteration count* — and it still does
  not move the exponent, because the per-application cost moves against
  it. Measured against the unmodified engine at three depths (value
  agreement over the full working precision in each case):

  | dps | baseline | accelerated | applications | cost/application | net |
  |---|---|---|---|---|---|
  | 200 | 20594 ms | 10140 ms | 2.16× fewer | 1.07× | **2.03×** |
  | 400 | 219766 ms | 99453 ms | 2.85× fewer | 1.29× | **2.21×** |
  | 800 | 2253812 ms | 9192078 ms | 2.06× fewer | 8.42× | **0.25×** |

  End to end between dps 200 and 400 the exponent falls only from 3.415
  to 3.294 — 0.12, not the 0.5 the iteration law alone would give. At
  dps 800 the gain is negative. Two things break there: the application
  saving stops growing, and the live Krylov basis (92 vectors × 4104
  coefficients × 811 digits ≈ 127 MB, against 2.5 MB at dps 200) has to
  survive every garbage collection. Capping the Krylov depth would bound
  that, but the superlinearity exists precisely because annihilated
  eigenvalues *stay* annihilated, so bounding depth bounds the gain.

  One implementation note worth keeping, because it halved the
  per-application cost: the matvec need not be the difference of two
  large quantities. Since the operator is affine, A·v = G(v) − G(0)
  (agrees with G(x₀+v) − G(x₀) to 7.06e−105), and G(0) is computed once
  per burst. G(0) must be evaluated *after* G(x₀), however — `staylor`
  caches the walk endpoints, and a walk taken from the zero polynomial
  locks in different branches, destroying the very branch-locking that
  makes the operator affine.
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
  gives the law above.

  That correction opened one last question, since a guard that could be bounded
  would change the exponent: is the loss per BLOCK? If chunking the points into
  N/K blocks capped the guard at g*K, the per-iteration cost would fall from
  Theta(p^2) to Theta(p log p) -- the Moroz mechanism without implementing
  Moroz. Measured (min digits of agreement between the tree and Horner at equal
  precision, over all 1280 points, guard 6*K in every cell):

  | deg P | K=64 | K=128 | K=256 | K=512 | K=1280 |
  |---|---|---|---|---|---|
  | 256 | 298.5 | exact | exact | exact | exact |
  | 512 | 272.8 | 286.9 | exact | exact | exact |
  | 1024 | 257.2 | 233.6 | 262.9 | exact | exact |
  | 2560 | 257.4 | 201.1 | 114.8 | 103.0 | exact |

  The threshold is `K >= deg(P)/2`, exactly. The guard requirement is set by the
  POLYNOMIAL DEGREE, not by the block size, so K cannot be capped below N/2 and
  the guard stays Theta(N). The only exact configuration measures 3x slower than
  Horner. The route is closed.
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

**Where the exponent actually stands.** The iterate carries N(p)*p ~ 5.3 p^2
digits, and the iteration count is Theta(p) -- measured, and closed from four
independent directions (extrapolation, Krylov, spectral, asymptotic form).
A caveat on the fourth of those and on the operator route generally: what was
measured is that the iteration matrix is not *banded* and has no global rank
deficit. Neither property is the one a fast direct solver needs -- a circulant
is not banded, and a hierarchically off-diagonal-low-rank (HSS/HODLR) matrix has
full global rank. Toeplitz structure is separately excluded, because the column
norms span 32 decades (0.0886 at mode 5 to 1.5e-33 at mode 300) where a Toeplitz
matrix would have near-equal ones. But **off-diagonal rank has not been
measured**, and it is the one structure class that would still collapse the
Theta(p) factor. That is an open question, not a closure. Any
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
