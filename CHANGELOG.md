# Changelog

## v0.1.1 — 2026-07-30

Three correctness fixes to the fork. Each returned a plausible-looking wrong
number with no error and no warning; all three are now either right or loud.

### Correctness fixes

- **Large towers no longer return a placeholder.** For base 2 above x ≈ 4.9,
  `sexp` returned `1.4426950408889634E400` — which is `1E400/ln 2`, an internal
  saturation sentinel — instead of `2^65536 ≈ 2.0035E19728`. `safefs()` also
  substitutes a `random()` imaginary part above 5E8, so repeated calls
  disagreed and one measured run returned a *negative* value, impossible for a
  real base-2 tetration. The reconstruction step inside `invabel` now uses the
  exact map: you get the value to full precision, or PARI raises `e_OVERFLOW`
  when the tower genuinely leaves the representable range. **Behaviour change** —
  code that appeared to work on huge arguments now raises. `safefs()` keeps its
  saturation inside `abel()`, where it is load-bearing.

- **Real bases 1 < b < e^(1/e) use a different construction now.** They went
  through the Kneser machinery and delivered roughly 15 correct digits no matter
  how many were requested (b = 1.2 at dps 60), and unusable output above
  dps ≈ 100 — silently. The unmodified original behaves the same way. These
  bases now use regular iteration at the real attracting fixed point
  (Koenigs–Schröder), the classical construction for this regime: 403 of 403
  requested digits, and dps 60 drops from 20.5 s to 15 ms.

  **New requirement:** in raw PARI/GP the base must be given *exactly* — write
  `sexpinit(6/5)`, not `sexpinit(1.2)`. The construction needs the base to
  roughly twice the target precision, and a decimal literal only carries the
  current `\p`; `sexpinit` refuses an inexact base instead of quietly computing
  with a slightly different one. The Python backend converts decimal *strings*;
  a Python `float` is refused for the same reason.

  In this regime only `sexp` and `slog` are defined. `abel`, `invabel`,
  `sexptaylor`, `slogtaylor` and `halfsexp` raise, because the Kneser
  quantities (`ct`, `circc`, `circr`, `k`, `lnb`, `rslog`, `L`, `L2`) still
  belong to the previously initialised base and would otherwise report *its*
  values without complaint.

- **Bases just above e^(1/e) no longer lose their digits.** At b = 1.4494 and
  dps 150 the engine printed 38 digits of which 1.837 were true — no error, no
  warning. The per-level degree truncation introduced in v0.1.0 measured its
  tolerance against the largest *coefficient* rather than the largest
  *contribution*; those coincide only while the coefficients decay, which fails
  in this band. Fixed: 147.410 true digits at dps 150 (0.98·dps), digit-for-digit
  equal to the unmodified original across all of the original's valid digits.

### Paper VI: anchor-pure base reconstruction

- **New paper** -- `papers/paper6_anchor_pure_reconstruction.pdf`, "Anchor-Pure
  Reconstruction of Regular Tetration: Riemann-Hilbert Selection, Analytic
  Base Flow, and Blind Validation". No Zenodo DOI yet; the series table says so
  rather than leaving the column looking archived.
- **New modules** `fatou_backend.reconstruction` (three modes: fast, validated,
  certified) and `fatou_backend.certification` (certificate auditor), with the
  console entry points `tetration-reconstruct` and `tetration-cert-audit`. The
  calculator never initialises the *target* base engine at run time -- that
  would make the result fast after caching, but no longer anchor-pure.
- **The `8r1` segment certificate is produced here, not imported**
  (`research/certification8r1/`). The trajectory is computed by this
  repository's own engine and matched the one shipped with the paper bit for
  bit -- `sha256(nodes) = da3898e5...` -- across a different fork revision, a
  different PARI binary and a different operating system. The large interval
  arrays (18.5 MB) are regenerable output and are deliberately not committed;
  the trajectory, the code and the certificate are.
- **A defect in the shipped certificate was fixed before anything was
  generated.** `python-flint` truncates every series operation at `ctx.cap`,
  default 10; none of the three producers set it, while their payload declared
  `"beta_series_order": 22`. Every beta model was a degree-9 polynomial.
  Coefficients carrying information went from 5700/12672 (45.0 %, exactly 10 of
  22) to 12540/12672 (99.0 %). `replay_8r1.py` now takes the order *from the
  payload*, so a declared/computed mismatch surfaces instead of passing.
  Consequence: the five interval values quoted in Paper VI come from the
  order-10 computation and differ from what this repository now produces.
- **That fix uncovered a second, deeper defect, and it is not repaired.** At the
  declared order the Hilbert producer stops with "stopped Koenigs recursion did
  not contract". Both endpoints of the trace integral are exact singularities
  of the recursion `w <- log_b(w)` (x = 0 gives y = 1, so w1 = 0; x = 1 gives
  y = b, so w2 = 0), the graded quadrature clusters nodes at both ends by
  design, and only the lower endpoint is special-cased. Measured against the
  closure's own admissibility criterion (Satz 1.3,
  `sum_{k>=1}|c_k| R^k < |c_0|` at R = 1), **55 of 5184 graded rows have no
  valid Taylor model of `log`** while the code adds the fixed constant
  `BETA_SERIES_REMAINDER = 2e-16` anyway. The ratios are identical at order 10
  and order 22, so the shipped certificate did not avoid this condition — it
  truncated below it. Raising `KOENIGS_MAX_DEPTH` would make the producer
  terminate (depth 495 suffices) without making the model valid, and is
  therefore not done. See `research/certification8r1/README.md`.
- **The graded Hilbert trace is now produced through the closure layer**
  (`certify_hilbert_derived_8r1.py`). `koenigs.koenigs_field` was written for
  this job -- its docstring says "used for the graded Hilbert trace" -- and was
  called by nothing. Every remainder is now derived; where the closure cannot
  prove one it writes `"u": null` with the reason instead of a number, and
  `replay_8r1.py` accepts both Hilbert schemas but stops before the quadrature
  when any row is unproved, so a `null` can never quietly become a zero. The
  producer shards by panel and node range (six concurrent shards, ~26 min
  instead of 2.5 h) and `merge_hilbert_derived_8r1.py` re-derives the global
  gates from the merged rows rather than trusting each shard's summary.
- **The endpoint limit is analytic, not numerical.** Γ's derived remainder
  degrades smoothly toward the upper endpoint -- 1.9e−16 at 1−x = 3e−14,
  4.8e−03 at 1e−17, construction failure below 5.2e−18 -- while `|y'|` holds at
  1.365, so the division is not the cause. A deliberately finer radius ladder
  buys only remainders of 1e−3 to 4, and recomputing at 448, 896 and 1792 bits
  leaves the admissibility ratio identical to five digits (2.2029e+02 at node
  862) while the ball radius falls to exactly 0. The state data carries 64
  significant digits, so it is not a data limit either.
- **Measured result of the full derived chain: 172 of 5184 graded rows could
  not be certified** -- 55 because no model exists on the panel at all (exactly
  the 55 the R = 1 criterion predicts) and 117 because no radius above 1 admits
  a Cauchy tail. All 172 sit at quadrature index 830-863, the top 34 of 864,
  and none at the lower end.
- **The upper endpoint model, derived** -- see
  `research/certification8r1/docs/UPPER_ENDPOINT_MODEL.md`. Gamma is the
  Koenigs-Abel velocity, it satisfies
  `Gamma(x) = c*T(x)*[Gamma(x-1) - T(x-1)]`, and `sexp_b(-1) = 0` for every
  base gives `Gamma(0) = 0` -- the lower endpoint model `certify_8r1.py`
  already asserts, here derived -- and one step further the value nobody wrote
  down: **`Gamma(1,beta) = -c*b = -e^(beta + e^beta)`**, exact and entire in
  beta, confirmed against the orbit to twelve digits. The relation needs only
  one logarithm, of an argument near `b`; the step that kills the orbit is the
  second one, whose argument goes to zero. Applied to the 288 rows with
  `1 - x <= 1e-16`, all 172 refusals disappear.
- **Cause of the asymmetry, corrected.** An earlier note here blamed
  `T(1,beta) = b(beta)` moving with beta. That reads the symptom right and the
  cause wrong: `sexp_beta(-1) = 0` for every beta, so the puncture is
  beta-stationary at both ends. The real cause is the state's degree-6 Lobatto
  interpolation error in beta, `eta = T(1,beta) - b(beta)`, which is 1e-64 at
  the panel edges (those *are* Lobatto nodes) and 1.5e-18 to 5.0e-18 at the
  centres. It displaces the zero of w2 to `1 - x* = eta/T'(1) = 1.105251e-18`,
  and quadrature node 862 sits at 1.113721e-18 -- the displaced zero lands on
  the node to three digits.
- **The segment certificate now closes at the declared order.**
  `replay_8r1.py` against the derived payload: PASS, Z1...Z7 PASS, `q < 1`,
  `Y + (q-1)r < 0`. The published values are about **eleven times tighter**
  than this chain can justify -- Y, r and p(r) all scale by the same factor
  10.890 (r follows from Y, and `p(r) = Y + (q-1)r`), while L and q do not move
  at all. The conclusion stands; the stated intervals do not reproduce. One leg
  is still on a set constant: the interval primitives come from
  `certify_8r1.py` with `BETA_SERIES_REMAINDER = 2e-16`, so the result is of
  mixed provenance.
- **Both legs of the `8r1` chain are now derived** (2026-09-06).
  `certify_derived_8r1.py` routes the interval primitives through the closure
  layer (`koenigs.koenigs_pair`, `tails.certified_tail`) with the same row
  layout as the constant producer, six panels in parallel; `replay_8r1.py`
  accepts both interval schemas, records the provenance of each leg in its
  report, and stops before the arithmetic when a payload carries unproved
  rows. The all-derived replay passes: `Y = 3.151553145348560e-9`,
  `L = 13613.30588114737`, `q = 0.1495393411945083`,
  `r = 4.63212659033378e-9`, `p(r) = -7.8788828633714e-10`, Z1...Z7 PASS --
  the same numbers as the mixed-provenance chain to eight digits, and still
  the factor 10.890 against Paper VI's order-10 table.
- **A soundness defect in the closure's tail rule, found and fixed.** With two
  interval producers that share no remainder logic, their enclosures were
  compared as values at five base samples per row: `tail_partial_T` disagreed
  on 254 of 2880 samples. The closure had charged the `d/dT` sensitivity tail
  with the value ratio `1/(c y_{j+1})`; the true ratio carries the factor
  `1 + c y_j theta` -- about twenty at the truncation index. At the worst row
  the first neglected term (8.29e-9) exceeded the claimed remainder (6.16e-10)
  by a factor of three. `tails.py` now uses `(1 + u_N)/(c e^{u_N})` for that
  sensitivity, with the tower's monotonicity above `log(2/c)` as a checked
  premise; `test_tails.py` computes the first neglected tower level explicitly
  and keeps the old rule as a negative control. The `d/dT'` sensitivity and the
  value were bounded correctly. After the fix all 34560 value comparisons
  between the two producers overlap. The closure document carries a dated
  erratum.
- **Z5 for derived payloads.** The replay's tail-convergence gate compared the
  last summed term with `TAIL_STOP = 1e-50`, the constant producer's stopping
  rule. The closure stops as soon as geometric domination is certified and
  folds the proved remainder into the arrays, so for the derived schema the
  gate now checks the certified ratio (`1 - max ratio_upper > 0`) and reports
  the rule it applied. The gate for the constant schema is unchanged.
- **Independent verification at the declared order.** The verdict file the
  README had promised now exists: `verify_primitives_8r1.py` regenerated the
  Koenigs and tail primitives at order 22 against the constant producer's
  arrays and passed; see `research/certification8r1/README.md` for the
  numbers.
- **The certificate template exists.**
  `research/certificates/paper_vi_certificate_template.json` is the schema-v2
  layout the auditor reads, written from the auditor rather than from any
  computed result: every number and hash is a `REQUIRED` placeholder, every
  claim flag is `false`, and a test binds its key tree to the auditor's own
  certificate fixture. `tetration-cert-audit` rejects it and lists what a
  producer must supply, which is what a template is for.
- Scope, stated plainly: this certifies **one** segment. The march segments,
  the segment-uniform gluing and the cut continuation are analytically
  specified and numerically open, and the closure document lists two of its own
  five points as unclosed. `fatou_backend.certification` refuses to report a
  final certified status by design.

### Verification

- The sub-eta test used to assert that the fork agrees with the unmodified
  original. That oracle was worthless: both implement the same construction, so
  they agree on its errors — which is precisely why the defect above survived
  v0.1.0. The fork test and the backend test now compare against an independent
  implementation (`research/tools/regular_subeta.py`).
- New test `test_near_eta_base_survives` pins the band just above e^(1/e).
- `sexp_e(0.5)` is now proven to **1973 digits** (dps 2000/2033 error-vector
  pair), up from 992.
- Tests: 96 passing, 4 skipped. `research/gate.py` still passes 11/11 with
  bit-identical values — the fixes touch regimes the gate does not cover and
  leave the ones it does cover untouched.

### Known limitations

- The truncation fix costs a measured **+9.7 %** on base e (median of five
  alternating runs per variant, ranges disjoint). It keeps more terms, which is
  the correction itself, not an implementation cost. Taken deliberately:
  correctness in a regime the old code destroyed outweighs 9.7 % on one base.
- Every cached fork state from v0.1.0 is invalid — the cache key hashes the
  engine file. Stale entries are unreachable rather than wrong, but they stay on
  disk until `state_cache.prune_cache(...)`.
- The vendored unmodified `fatou.gp` is third-party code and still carries all
  three defects. Use the fork.
- Sub-eta support covers `sexp` and `slog`, not the Taylor/Abel entry points.

### Corrections to v0.1.0 claims

- v0.1.0 listed the sub-eta regime under end-to-end coverage "against the
  unmodified original". That verification is retracted; see above.
- The Θ(p³) floor in v0.1.0's *Known limitations* holds for any scheme that
  touches its state once per iteration. Off-diagonal rank has never been
  measured, and it is the one structure class that could still collapse the
  Θ(p) factor — an open question, not a closed one (`research/METHODS.md` §3).

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
  tiers and `sexp_2(0.5)` to 497 (error vector **and** engine diversity). A
  dps 1200/1265 pair extending this to ~1165 digits exists as a candidate
  (`research/tools/error_vector_pair.py`) and is not yet part of the frozen
  ladder.
- `research/tools/digits_vs_reference.py` measures true digits against those
  references and refuses to report agreement past the proven ceiling. Its
  `--pin` option pins the run to a fixed core set and verifies the mask took
  effect, refusing to report timings it could not actually pin.
- The DFT/twiddle layer is tested directly (`tests/test_transform_layer.py`):
  against a naive O(n²) DFT, against high-precision truth for the power
  tables, and by a source check that no `powers()` call survives at a twiddle
  site — the construction that carried ~n ulp instead of ~log₂(n).
- Tests: 83 passing, 4 skipped (`python -m pytest tests/`), including
  end-to-end coverage of the fork itself (engine diversity against the
  unmodified original, the sub-eta regime, and the ≥0.95·dps calibration
  claim), the two console entry points, and the DFT/twiddle layer.
- CI (`.github/workflows/tests.yml`) runs the fast suite on Linux, where
  PARI/GP is `gp` on PATH — the configuration most likely to break and
  least likely to be noticed on the author's machine — and asserts that a
  built wheel actually contains its runtime data.

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

### Interface notes

- `FatouGP` and the CLI return values at full precision regardless of the
  caller's `mp.mp.dps`; that global still governs how they *print*, so use
  `mp.nstr(v, n)` or raise it to see the digits.
- Unsupported bases (b ≤ 0, 0 < b < 1, b = 1) are rejected immediately with an
  explanatory error instead of hanging until `init_timeout`.
- `eval_batch`'s `digits=` parameter is gone: it was accepted and ignored.
- The state cache has no automatic eviction; `state_cache.prune_cache(...)`
  and `state_cache.cache_size()` are there when the directory grows (measured:
  159 states = 72 MB).

### Known limitations

- **The depth-scaling exponent is ~4.0–4.16 and this release does not change
  it.** All optimizations here are constant-factor by construction; the local
  exponents confirm it (300→400: 3.31 before, 3.24 after). The structural reason
  is in `research/METHODS.md` §3: with Θ(p) Picard iterations each touching
  Θ(p²) digits, Θ(p³) is a floor for this discretization.
- Measured core-pinned and bracketed, on a session whose bracket showed 0.991
  drift (i.e. the machine held still): **4.159** for 520→1200, 4.151 for
  520→1265, 4.025 for 1200→1265. An earlier bracketed attempt gave 4.35 but ran
  on hardware sitting 1.4–1.7× off its own best state; it is superseded, since
  the exponent grows with depth and 520→1200 cannot be *below* 520→1020.
  Net: unchanged from the historical 4.13.
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
