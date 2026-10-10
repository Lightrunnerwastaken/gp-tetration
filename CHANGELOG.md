# Changelog

## Unreleased

### Pentation at the engine's precision, for every base the engine runs

- **`pentinit` / `pent` / `ipent` now carry the engine's digits.** The
  original builds the Poincare function of slog at the negative fixed point of
  sexp with a fixed 21-term series (`formalischroder(tslog,21)`), so its
  pentation delivered ~30 digits regardless of the requested precision -- at
  `\p 60`, `pent(0)-1 = 6.5e-31`. The fork now grows the series with the
  working precision (Newton doubling instead of the original's O(N^4)
  term-by-term composition) and chooses the slog sampling from the measured
  decay of the series. At dps 40 the defining identities hold to 1e-48 ..
  1e-57; dps 60 against dps 100 agree to 1.1e-65 relative on `pent(1/2)`.
- **New bases.** The construction now goes through `slog()` instead of
  `slogtaylor()`, so it also runs for real 1 < b < e^(1/e) (regular tetration,
  where `slogtaylor` is guarded) -- the original cannot run there at all -- and
  for complex bases, where the original's error estimate (`z1<1`) failed.
  Checked: 21/20, 6/5, 7/5, 36/25, 29/20, 3/2, 17/10, 2, e, 10, 100, 2+i/2.
- **A precision trap, fixed before it shipped.** For base 100 (lambda = 1100,
  fixed point 2e-4 above the singularity at -2) the Newton solve leaves the
  highest Poincare coefficients with 19 of 57 digits. Their contribution at
  the evaluation radius is below 1e-65, but PARI carries precision through
  Horner relatively, and `pent(0)-1` came out 3.6e-21. Terms whose error
  contribution is provably below the target are now padded (as `regsigma`
  already does); otherwise `pentinit` raises.
- **`pentinit` reports a real error estimate.** The original returned 0
  whenever `sexpinit` was skipped. Now: `pent(0)=1`, an inverse roundtrip and
  the same values over two reduction depths.
- Unchanged: `sexp`/`slog` are bit-identical to the previous fork (checked at
  dps 60 on e, 1.2 and 2+i/2), the original `fatou.gp` keeps its version, and
  Sheldon's hexation (`hexinit`, `hex`, `ihex`) runs on top of the new
  pentation as before (base 2: `hex(1.5) = 2.67289396602652...`, its own
  estimate ~8 digits).
- New tests: `tests/test_pentation.py` -- identities on both sides of eta, a
  large and a complex base; dps 40 against dps 60; Fox (2007) to 1e-20;
  Levenstein (2020) for base 2; the original's pentation to 1e-25.

### Hexation at full precision for real 1 < b < eta5

- **`hexinit` now takes the regular route where one exists.** For real bases
  below eta5 = 1.7373559... the pentation crosses y = x twice and the lower
  crossing a is attracting (b = 1.7: a = 2.3876372817, P'(a) = 0.449). There
  the hexation is the regular iteration at a -- the construction the fork
  already uses for tetration below eta: `hex(z) = Psi_H(kappa*mu^z)`,
  `hex(0) = 1`, `hex(z+1) = pent(hex(z))`, `hex(z) -> a`. At dps 40 the
  identities hold to 1e-47 .. 1e-57 (1.2, 1.41, 1.5, 1.7, 1.73), and dps 40
  against dps 60 agree to 1e-46 .. 1e-57.
- **Above eta5** (2, e, 10) the pentation has no real fixed point; see the
  Kneser hexation below.
- `hexinit(B, n, skip)` gains `skip` (reuse the initialised engine, as
  `pentinit`); `pentinit` resets the hexation state of a previous base.
- New tests: `tests/test_hexation.py`.
- Found on the way: GP's `local()` is dynamically scoped. A helper that
  declared `local(L)` handed its own `L` to every engine function it called,
  and the Kneser engine reads the global `L` (its fixed point) -- slog was
  silently wrong for every base above eta, right below eta (regular path, no
  `L`). The new pentation/hexation code uses `my()` throughout.

### Kneser hexation for real bases above eta5

- **`hexinit` now builds the Kneser hexation for b > eta5** instead of Sheldon
  Levenstein's ~8-digit matrix approximation, which stays as the fallback.
  At the complex fixed point p of the pentation (b = 2: p = -1.99295 +
  1.44513 i) the Abel function of the hexation is `alpha_p + theta(alpha_p)`
  with a 1-periodic, decaying theta; between p and the real axis it is a real
  Taylor series around Re p. One least-squares problem couples both: theta
  matching near p, the functional equation near the real axis (every arc
  point gets an equation), the normalisation alpha(1) = 0. Solved by modified
  Gram-Schmidt -- PARI's `matqr` only takes square matrices, and the normal
  equations were singular for base e.
- **The digit count is measured, not assumed.** Every solve is repeated on a
  second arc with its own points and scaling; `hexerr` is their disagreement.
  A small residual alone proves nothing here: with arc points left without an
  equation the system was not unique, and a series 2e-21 off satisfied every
  equation to 1e-32.
- Measured: base 2, `hex(1.5) = 2.6728939593426194891907964716798365760`;
  `\p 20`: 32.6 digits in 14 s, dps 45: 43.4 digits in 51 s, the two agree to
  2e-32. Sheldon's value is 7e-9 off (his own estimate: 8.1 digits), Levenstein
  2020 quotes 2.6729. Base e, `hex(1.5) = 9.645503939964672449631960...`: about
  26 digits, at dps 20 and at dps 45 alike (residual 2e-26 at dps 45, cause
  open; the error vector reports it).
- **Base 10 was slow (superseded by exp-081 below).** p sits close to the
  axis (Im p = 0.87), the arc has to reach 0.85 Im p, K = 740: 4.3 min at
  `\p 20` for 13.5 digits. `hexinit` skips Kneser above `knmaxk = 600` terms;
  set `knmaxk = 900` to force it.
- Real bases only. Complex arguments of `hex`/`ihex` work near the real axis
  (Newton inside the fitted disk).
- Tests: `tests/test_hexation.py::KneserHexationTests` (base 2 by default,
  base e with `FATOU_BACKEND_RUN_SLOW=1`).

### Kneser hexation: centre planned, radius measured (exp-081)

- **The expansion centre is no longer fixed at Re p.** For base 10, p lies
  close to the axis (Im p = 0.87); a disc around Re p had to reach 0.85 Im p to
  hold one period of the real axis with its images (K = 740, 4.3 min at `\p 20`,
  13.5 digits). `knhexinit` now plans centre and radius together, in pairs of
  neighbouring radii (the error vector needs the second arc).
- **Images must stay inside the arc** (reach <= 1.02 rad, one relaxed pass at
  1.15 if nothing qualifies). Images outside forced A to extrapolate: that was
  the old base-e ceiling of ~26 digits (found by the tetration-graph session).
- **K is replanned from the measured decay** of the Taylor coefficients after
  the first solve -- the planner's radius |c - p| is too optimistic; alpha has
  nearer singularities (base e near 0.33 +- 1.37 i).
- Robustness: `ipent` stops after 60 slog steps (it looped forever for some
  complex points and hung the planning), every planning candidate runs under
  `alarm(1)`, and slow candidates are skipped; an inner `iferr` no longer
  swallows the alarm.

| base | dps | before (digits, time) | now (digits, time) |
|---|---|---|---|
| 2 | 20 | 32.6, 14 s | 31.9, 10 s |
| 2 | 45 | 43.4, 51 s | 47.6, 27 s |
| e | 20 | 29.6, 25 s | 30.8, 42 s |
| e | 45 | 25.7, 74 s | 42.2, 42 s |
| 10 | 20 | 13.5, 258 s | 33.4, 21 s |
| 10 | 45 | -- | 43.6, 70 s |
| 50 | 20 | Sheldon, 0.2 | 21.5, 64 s |

  Base 100 still fails: the search for the complex fixed point of the
  pentation (ipent iteration from 0.5) does not converge there.

### Kneser hexation: fewer theta harmonics, profiled (exp-082)

- **Fewer harmonics, faster and more accurate.** The N formula assumed theta
  decays at height 0.6*(htop - hmin); the corners of the real curve lie
  below hmin and the c_n fall on top of that, so the effective height is
  ~0.3 larger. Too many harmonics hurt (base e, dps 45, same arc: N = 80
  residual 9e-24, 49 -> 2e-42, 35 -> 7e-47, 25 -> 8e-47).
- K is replanned from a fit log|a_k| ~ a + b k + c log k instead of a
  two-point slope.
- Profile (base e, dps 45, 58 s before): the least-squares solves were 47 s.
  Tried and rejected: classical Gram-Schmidt via matrix-vector products
  (slower -- the products also run over the empty columns), one
  orthogonalisation pass (-10 digits), 1.5x instead of 2x working precision
  (no gain: PARI computes in whole words).
- Now: dps 45 -- base 2 46.3 digits in 25 s, e 44.4 in 33 s, 10 43.3 in 63 s;
  dps 20 -- 2 31.8 in 8 s, e 30.7 in 31 s, 10 34.4 in 20 s.
- **The ceiling is the input, not the solve.** At dps 45 the residual stays at
  ~1e-47 for larger K and for 3x working precision. The pentation values the
  equations are built from are accurate to 1e-51..1e-52, with spikes to
  4e-49 that come from the engine's sexp (1e-52..1e-55, spikes at particular
  fractional parts of x). More hexation digits need a more accurate sexp.

### Kneser hexation: wider planning for large bases (exp-087)

- The relaxed planning pass -- reached only by bases where no arc keeps its
  images inside, b >~ 20 -- also tries radii 1.45 and 1.6 Im p and pairs equal
  radii around neighbouring centres; the second arc is solved around its own
  centre. Base 100, dps 20, K = 300: centre Re p + 1.5 Im p, radius 1.45 gives
  residual 3e-31, the old main arc (Re p + Im p, radius 1.3) 1e-13 -- it was
  outside the radius list and had no neighbouring radius.
- The shifted sibling arc (exp-085) now runs before the own-K / half-K retry
  of the second arc: for b = 100 at dps 20, K = 600 on the second arc cost
  64 s for 7e-29, and the sibling afterwards gave 2e-32.

| base | dps | exp-086 (digits, time) | now |
|---|---|---|---|
| 100 | 20 | 22.2, 98 s | 31.6, 119 s |
| 100 | 16 | 12.1, 15 s | 12.7, 17 s |
| 50 | 20 | 32.7, 108 s | 32.6, 60 s |
| 10 | 45 | 52.9, 86 s | 53.0, 83 s |

  Unchanged: b = 2 and e at dps 45 (51.0 in 35 s, 53.2 in 53 s), b = 28.06 at dps 16 (13.0 in 18 s).
- Tests: full suite 123 passed / 6 skipped; level tests incl. slow 27 passed.

### Kneser hexation: the second arc gets its own K (exp-086)

- The half-K retry of exp-084 is right only when the second arc's images lie
  well outside it (its error grows like (rho/rad)^K). When they lie almost
  inside (rho/rad <= 1.05), the arc is short of terms: it now gets its own K
  from its own coefficient decay (`knfit`, capped at knmaxk), and if it then
  beats the main solution the two swap roles.
- Base 50, dps 20 -- arc 1.3 Im p (rho/rad 1.015): K 221 residual 6e-23,
  K 330 1e-32; arc 1.15 (rho/rad 1.18): K 221 6e-32, K 330 5e-31. No
  cancellation (sum |a_k| = 4.2). Now arc 1.3 with K = 494 (2e-33) is the
  main solution and arc 1.15 checks it: **23.1 -> 32.7 digits**, 108 s.
- Unchanged: b = 10 at dps 45 (52.9), b = 100 at dps 20 (22.2 -- there the
  main arc itself would need K > knmaxk), b = 2 and e at dps 45 (51.0 in 27 s, 53.2 in 51 s).
- Tests: full suite 123 passed / 6 skipped; level tests incl. slow 27 passed.

### Kneser hexation: a shifted sibling arc certifies base 10 (exp-085)

- If the second arc stays much worse than the main one even after the half-K
  retry (images outside, exp-084), one more solution is computed on an arc of
  the main radius around a centre shifted by +0.125, +0.25 or -0.125 Im p (at
  most three tries) -- an independent solution whose images stay inside.
  The error vector now evaluates each solution around its own centre.
- Base 10, dps 45: **46.5 -> 52.9 certified digits** (|s1 - s2| = 2.3e-54),
  86 s. Shifted centres measured: +0.125 Im p residual 1e-53, +0.25 1e-53,
  -0.125 2e-44, -0.25 1e-34.
- Where it does not trigger nothing changes: b = 2 and e at dps 45, b = 100
  (both arcs within a factor 100), b = 50 at dps 20 (relaxed arcs, no
  shifted plan passes the reach test).
- Tests: full suite 123 passed / 6 skipped; level tests incl. slow 27 passed.

### Kneser hexation: base 100, the better arc leads (exp-084)

- **Base 100 works.** The series of the pentation around its complex fixed
  point p has radius 0.0075 there (b = 50: 0.013), but the sampling radius
  was halved only down to 1/128. Now down to 1/2048 (b = 100 needs 1/256).
  The ipent iteration did find p = -2.6587 + 0.4426i; the earlier note
  "fails in the fixed-point search" named the wrong step.
- The plan may ask for up to 2 knmaxk: it overestimates K (b = 28: plan 575,
  287 sufficed), and the pre-solve with K0 = plan/2 <= knmaxk measures what
  is needed (b = 100, dps 20: plan 1094).
- **The better arc is the main solution.** The error vector took the
  residual of the first arc even when its neighbour was much better (b = 10,
  dps 45: 1.1e-46 against 1.2e-53).
- **Images outside an arc amplify its error** like (rho/rad)^K: b = 10,
  dps 45, arc 1.0 Im p has rho/rad = 1.06 and residual 1e-46 at K = 352,
  7e-44 at K = 460; arc 1.15 keeps its images inside and gives 1e-53 for every
  K and N tried. The planner's 30-point reach misses this (1.0014). If the
  second arc is much worse, it is re-solved with half K and the better one
  kept (b = 100, dps 20: 3e-9 -> 4e-22). A third arc further out does not
  help -- it needs more K (b = 10, arc 1.32, K = 352: 4e-32).

| base | dps | exp-083 (digits, time) | now |
|---|---|---|---|
| 100 | 16 | fails | 12.1, 15 s |
| 100 | 20 | fails | 22.2, 100 s |
| 10 | 45 | 46.0, 64-82 s | 46.5, 48 s |
| 50 | 20 | 23.1, 71 s | 23.1, 72 s |
| 2 | 45 | 51.0, 30 s | 51.0, 21 s |
| e | 45 | 53.2, 65 s | 53.2, 47 s |

  Base 10 at dps 45 stays at ~46.5 certified digits: the error vector needs a
  second solution of similar quality, and on the smaller arc the images leave
  the disc. The main solution has residual 1e-53 and round trip 1e-56, so it
  is probably better than certified.
- Calculator: hexation also for b = 100 (11.6 digits) -- the whole base
  grid 1.5 ... 100 now has both levels.
- Tests: full suite 123 passed / 6 skipped; level tests incl. slow 27 passed.

### sexp spikes removed; Kneser K from a pre-solve; hexation tables to b = 50 (exp-083)

- **sexp spikes (engine).** `betterest` stopped refining at
  10^(-precis/1.2), so at particular fractional parts of x sexp held only
  ~5/6 of the working precision (dps 60: worst point 2.5e-64, median
  ~4e-74). It now stops at 10^(-precis+2): worst 1.3e-73, median unchanged,
  ~25 % more time per sexp value. Gate PASS (6 min 22 s).
- Pentation, base e, dps 45 against dps 70: error 4e-52 ... 5e-49 ->
  4e-55 ... 1e-52; ipent(pent(x)) - x 7e-52 ... 4e-49 -> 4e-57 ... 3e-56.
  The Kneser hexation was bounded by exactly these inputs (exp-082).
- **K from a pre-solve.** A solve with half the planned K (~1/8 of the cost)
  measures the coefficient decay (`knfit`); K = 1.1 x what the fit needs for
  target + 5 digits, where before the first full solve was often redone
  with a larger K. If the fit sees no decay, the coefficients already sit at
  their floor and the pre-solve is kept (b = 21.7 ... 37.1 at dps 16: plan
  K = 575, 90 s -> 12.5 s, same digits). If the larger K makes the solve
  worse, the pre-solve is kept as well: arc points lie at |u| = rad/rho < 1,
  images up to |u| = 1, conditioning ~ (rho/rad)^K (b = 50, dps 20: K = 600
  residual 1e-13, K0 = 221 6e-32). A fit K above `knmaxk` is capped instead
  of falling back to Sheldon's approximation (0.2 digits at b = 50).

| base | dps | exp-082 (digits, time) | now |
|---|---|---|---|
| 2 | 45 | 46.3, 25 s | 51.0, 30 s |
| e | 45 | 44.4, 33 s | 53.2, 65 s |
| 10 | 45 | 43.3, 63 s | 46.0, 64-82 s |
| 2 | 20 | 31.8, 8 s | 32.1, 8 s |
| e | 20 | 30.7, 31 s | 34.1, 24 s |
| 10 | 20 | 34.4, 20 s | 34.5, 13 s |
| 50 | 20 | 21.5, 64 s (exp-081) | 23.1, 71 s |

  Times measured with other jobs on the machine (+-30 %).
- `pentx` / `ipentx` / `hexx` / `ihexx` raise an error beyond the overflow
  instead of returning a zero without valid digits (b = 100, dps 16:
  pent(1.96875) = 0.E495855885102019914, which Python read as 0).
- **Calculator: hexation up to b = 50.3** (7 new bases 11.2 ... 50.3,
  12.6-14.5 digits each). Along x the tables are now interpolated by a
  monotone cubic (Catmull-Rom tangents limited after Fritsch-Carlson):
  even in symlog space the curves explode (b = 21.66: 4.5e7, 7.8e20, 5.9e220
  one step apart) and Catmull-Rom undershot -- hex(1.2) came out as 1875
  instead of ~9e8. Tables end at their first overflow; a stray 0 behind the
  gap stood in 14 tables (10 pentation, 4 hexation).
- Tests: `test_overflow_is_an_error_not_zero`; full suite 123 passed / 6 skipped,
  level tests incl. slow 27 passed.

### Pentation and hexation from Python and the CLI

- `FatouGP.pent / ipent / hex / ihex` (and `_batch` forms, also on sessions),
  `pentation_digits(base)` and `hexation_digits(base)`; CLI commands `pent`,
  `ipent`, `hex`, `ihex` (`fatou-cli pent --base 2 --values 1.5 --fatou-gp fork`).
  Fork only -- the original is refused with a message, since its pentation
  carries ~30 digits at any precision.
- The fork's `pentx/ipentx/hexx/ihexx` initialise the level on first use per
  (base, precision) inside the GP process (`pentensure`, `hexensure`), so each
  expression stands on its own and a pool of workers works unchanged.
  `hexinit` is split into `pentinit` + `hexsetup`.
- The first hexation call above eta5 builds the Kneser solution (~15-25 s at
  dps 20, ~1 min at dps 45).
- Calculator: `tetration_calculator.html` gains a level switch (tetration /
  pentation / hexation). The curves are tables computed with the fork by
  `research/tools/calc_export_levels.py` (dps 16, step 1/32): pentation for
  77 bases 1.5 ... 100, hexation for 49 bases 1.5 ... 100 (incl. 2, e, 3, 10;
  12-15 digits each). Monotone cubic in symlog space along x (exp-083),
  linear in lnln b between bases, never across eta5, no extrapolation. At
  grid bases the page
  matches the fork to ~1e-13: 2^^^1.5 = 2.70289057146015, 2^^^^1.5 =
  2.67289395934264, e^^^^1.5 = 9.645503939965.
- Tests: `tests/test_levels_api.py`.

### State cache: restored engines are bit-identical to fresh ones

- **A cache-restored engine returned ~15 fewer digits.** At dps 40, base e,
  against a fresh dps-100 run: fresh 54-56 correct digits, restored 40-41
  (dps 60: 73-75 against 60-61). The internal values were bit-identical
  (checked with `binary()`); only the printing differed. `default(realprecision,
  n)` also sets how many digits `print()` emits, the worker reads results back
  from that text, and the fork's `loop()` ends with `default(realprecision,
  precis)` -- 57 digits at dps 40. The restore path skips `sexpinit` and kept
  the init line's 40, and `\uv` lists variables, not GP defaults. The cache
  now stores `realbitprecision`, `seriesprecision` and `format` in the JSON
  sidecar and restores them. Entries written before this fix lack them and are
  rebuilt once.
- **The cache never hit at mpmath's default precision.** The anchor check
  parsed the stored anchor at the caller's `mp.mp.dps`; at the default 15 it
  rounded to a double, missed the 1e-35 tolerance, and every entry was silently
  rebuilt. The cache only ever worked for callers who had raised `mp.mp.dps`
  -- which is also why results depended on it. The comparison now runs at the
  stored anchor's own precision.
- Checked fresh against restored, bit for bit: fork at dps 40 (e, 1.2, 2,
  2+i/2) and dps 60 (e), original at dps 40 (e). New test:
  `tests/test_state_cache.py::CachedEqualsFreshTests` (e and 1.2, at mpmath's
  default precision).

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
