# Paper VI — the `8r1` segment certificate

This directory produces the interval certificate for segment `8r1` of the
anchor-pure base reconstruction (Paper VI,
[`papers/paper6_anchor_pure_reconstruction.pdf`](../../papers/paper6_anchor_pure_reconstruction.pdf)).

The segment covers β ∈ [7β₀/8, 11β₀/12] with β₀ = log log 2, that is bases
b ≈ 2.043 … 2.060 — comfortably inside the Kneser regime, far from
η = e^(1/e).

## The data here is produced, not imported

The trajectory is computed **by this repository's own engine**. That is the
whole point: a certificate you cannot regenerate is a claim, not a
certificate.

```bash
cd research/certification8r1

# 1. trajectory from our fork          (~17 min)
python produce_8r1.py \
    --output trajectory_8r1.json \
    --gp-exe "C:/Program Files (x86)/Pari64-2-17-3/gp.exe" \
    --fatou-gp fork

# 2. interval primitives -- derived remainders, six panels in parallel
#    (~14 min wall; this is the proof path)
python certify_derived_8r1.py \
    --trajectory trajectory_8r1.json \
    --output segment_8r1_derived_interval_primitives.json --jobs 6

#    the constant-remainder producer (~7 min), kept because two producers
#    that share no remainder logic are the cross-check that found the
#    tail defect described below
python certify_8r1.py \
    --trajectory trajectory_8r1.json \
    --output segment_8r1_interval_primitives.json

# 3. Hilbert/trace primitives -- derived remainders, six shards in parallel
for p in 0 1 2 3 4 5; do
  python certify_hilbert_derived_8r1.py \
      --trajectory trajectory_8r1.json \
      --output "shards/panel$p.json" \
      --panel-start "$p" --panels "$((p+1))" &
done; wait                                     # ~23 min wall, not 2.5 h

python merge_hilbert_derived_8r1.py \
    $(for p in 0 1 2 3 4 5; do printf -- '--shard shards/panel%s.json ' "$p"; done) \
    --output segment_8r1_hilbert_derived.json

# the constant-remainder producer, kept for comparison; it does not run
# to completion at the declared order -- see below
# python certify_hilbert_8r1.py \
#     --trajectory trajectory_8r1.json \
#     --output segment_8r1_hilbert_primitives.json

# 4. independent replay of the stored arrays -- both legs derived
python replay_8r1.py \
    --primitives segment_8r1_derived_interval_primitives.json \
    --hilbert    segment_8r1_hilbert_derived.json \
    --trajectory trajectory_8r1.json \
    --report     segment_8r1_replay_all_derived_report.json

#    the same replay with the constant-remainder interval leg (mixed
#    provenance); its report is kept as the cross-check it is
python replay_8r1.py \
    --primitives segment_8r1_interval_primitives.json \
    --hilbert    segment_8r1_hilbert_derived.json \
    --trajectory trajectory_8r1.json \
    --report     segment_8r1_replay_derived_report.json

# 5. independent regeneration of the Koenigs and tail primitives (~15 min)
python verify_primitives_8r1.py \
    --trajectory        trajectory_8r1.json \
    --stored-primitives segment_8r1_interval_primitives.json \
    --output            segment_8r1_derived_primitives.json \
    --report            segment_8r1_verification_report.json
```

Step 4 **passes**, at the declared order 22, with every remainder derived by
the closure layer on both legs and the derived upper endpoint model of
[`docs/UPPER_ENDPOINT_MODEL.md`](docs/UPPER_ENDPOINT_MODEL.md)
(`segment_8r1_replay_all_derived_report.json`, 2026-09-06):

```
8r1 array replay: PASS
Y    = [ 3.151553145348560e-9 +/- 8.67e-25]
L    = [13613.30588114737     +/- 8.10e-12]
q    = [ 0.1495393411945083   +/- 7.33e-17]
r    = [ 4.63212659033378e-9  +/- 4.40e-24]
p(r) = [-7.8788828633714e-10  +/- 1.86e-24]
Z1 … Z7: PASS
```

with `q < 1` and `Y + (q−1)·r = −7.8788828633714e−10 < 0`. Without the
endpoint model the same chain rejects, and should:

```
issue: Hilbert payload carries 172 unproved rows;
       the graded Hilbert term is undefined and the replay stops
```

**Both legs of this chain now go through the closure layer.** Until
2026-09-06 the interval primitives came from `certify_8r1.py`, which charges
`BETA_SERIES_REMAINDER = 2e-16` at every node, and only the Hilbert leg was
derived — a verdict of *mixed provenance*. `certify_derived_8r1.py` routes the
interval leg through `koenigs.koenigs_pair` and `tails.certified_tail` with the
same row layout, and the replay report records which producer each leg came
from. Running the two interval producers side by side is also what exposed a
defect in the closure's tail rule (§ *A defect in the closure's tail rule*
below); the numbers above are from the chain with that defect fixed.

Of what steps 2–5 write, only the large primitive arrays (13–20 MB each) are
excluded by `.gitignore`; they are regenerable in minutes. The *reports* are
committed, because they are the certificate: the replay reports carry Y, L,
q, r, p(r) and Z1…Z7, and the verification report carries the independent
PASS. A repository with the code but not the result would be a recipe, not a
certificate. Requires `pip install python-flint` (the `cert` extra).

## Reproduction result, 2026-07-31

The trajectory produced here matched the one shipped with the paper **bit for
bit**:

```
sha256(nodes) = da3898e549bf88743df5a1295093626e09db197a7addca3a87970269540ed791
```

on both sides — across a *different fork revision* (ours carries exp-071 and
exp-074, 116218 bytes against 101793), a *different PARI binary*, and a
*different operating system*. The provenance manifest records those
differences rather than hiding them, which is why the payload hashes differ
while the numbers do not. That is the useful kind of agreement: not a replay
of the same environment, but the same result from a different one.

## A defect in the shipped certificate, and its fix

`python-flint` truncates **every** series operation at the global length
`ctx.cap`, whose default is `10`. None of `certify_8r1.py`,
`certify_hilbert_8r1.py` or `replay_8r1.py` ever set it, while their payload
declared `"beta_series_order": 22`. Every β model was therefore a degree-9
polynomial. This was found by the closure work
([`docs/PAPER_VI_8R1_CLOSURE.md`](docs/PAPER_VI_8R1_CLOSURE.md)) and confirmed
independently here (flint 0.9.0, `ctx.cap == 10` on this machine).

The producers now set `ctx.cap`; `replay_8r1.py` takes it *from the payload*
rather than from a constant, so a mismatch between declared and computed order
surfaces instead of passing unnoticed.

An earlier revision of this file claimed that fix covered all three producers.
It did not: `certify_hilbert_8r1.py` was still missing `ctx.cap` entirely, and
`certify_8r1.py` sets it inside `main()`, so importing that module does not
execute the assignment. The Hilbert primitives therefore stayed at order 10
while their payload declared 22 — invisible, because `encode_series` pads
shorter series with `acb(0)` up to 22 and the replay only checks the *length*.
It reaches the certificate through `hilbert[] → selected → G → residual → Y`.
Fixed 2026-08-03; every number in this file postdates that fix.

Measured effect on the Γ arrays of the interval primitives:

| | Γ coefficients carrying information |
|---|---|
| as shipped | 5700 / 12672 — **45.0 %** |
| after the fix | 12540 / 12672 — **99.0 %** |

45 % is exactly 10 of 22. The coefficients continue smoothly past the old
wall — Γ[9] = −2.52e−25 (bit-identical in both runs), Γ[10] = −3.79e−28,
Γ[11] = 2.31e−31, Γ[21] = 1.55e−57, a geometric decay of about three orders
per index. The twelve missing coefficients were real information, not
negligible padding.

## What the fix uncovered: the graded endpoint is singular

At the declared order the Hilbert producer **does not run**. It stops with

```
RuntimeError: stopped Koenigs recursion did not contract
```

in panel 1, and that failure is the honest one. Both endpoints of the trace
integral are exact singularities of the Koenigs recursion `w ← log_b(w)`:

| endpoint | state | consequence |
|---|---|---|
| x = 0 | y = 1 | w₁ = log_b(1) = 0 |
| x = 1 | y = b | w₁ = 1, w₂ = log_b(1) = 0 |

and the next step takes `log` of that zero. `certify_8r1.py` special-cases the
lower endpoint (`x_index == 0` → `endpoint_stopped_model`); nothing handles the
upper one. The graded Gauss–Legendre grid clusters nodes at both ends — that is
what grading is for — so the last nodes of every panel sit arbitrarily close to
a branch point.

The closure states the right criterion itself. `PAPER_VI_8R1_CLOSURE.md`
Satz 1.3 admits a Taylor model of `log P` on the disc of radius R only when

```
    Σ_{k≥1} |c_k| R^k  <  |c_0|          (which proves P ≠ 0 on D̄_R)
```

Measured against the β model of w₂ at R = 1, i.e. over the panel, that
criterion fails on 55 of the 5184 rows — worst ratio 220, at depth 3
throughout. Running the full chain through the closure layer (below) refuses
**172 of 5184**, and the split confirms the diagnosis exactly:

| refusal | rows | meaning |
|---|---|---|
| `polynomial part is not bounded away from zero on the disc` | 55 | no model on the panel itself; the same 55 the R = 1 test predicts |
| `no admissible Cauchy radius for the tail` | 117 | no R > 1, so no Cauchy tail to bound the remainder |
| **total** | **172 of 5184** | 23 / 32 / 34 / 33 / 33 / 17 per panel |

Every one of them lies at quadrature index 830–863, the top 34 of 864. **None
at the lower end.**

`certify_8r1.py` nevertheless adds the fixed constant
`BETA_SERIES_REMAINDER = 2e-16` at every node. Where no admissible R ≥ 1
exists there is no Cauchy tail to bound, and a constant is not an enclosure.

Two points matter for how this is read:

* **The order-10 run had the same defect.** The measured ratios above are
  identical to all printed digits at order 10 and at order 22 — this is a
  property of the data and the scheme, not of the truncation. Order 10 merely
  kept the coefficient blow-up at 10¹³ instead of 10³², small enough that the
  terminal-distance test still slipped under `KOENIGS_STOP = 1e-16` within
  `KOENIGS_MAX_DEPTH = 360`. The shipped certificate did not detect the
  condition; it truncated below it.
* **Raising the depth limit would hide it again.** At node 862 the box does
  reach 1e−16, at depth 495. That would make the producer terminate without
  making the model valid, so it is the wrong fix.

These rows reach the certificate through
`hilbert[] → selected → G → residual → Y`.

### Routed through the closure layer

`koenigs.koenigs_field` exists in the closure for exactly this job — its
docstring says "used for the graded Hilbert trace" — and was called by nothing.
`certify_hilbert_derived_8r1.py` now calls it: every remainder is derived, and
where the closure cannot prove one it refuses instead of emitting a number.
A refused row is written as `"u": null` with its reason, and `replay_8r1.py`
stops before the quadrature rather than letting a `null` become a zero — the
graded Hilbert sum needs every row, so a missing one is an undefined term, not
a small gap to widen.

The producer shards by panel and node range, so the chain runs as six
concurrent panel shards (~26 min) instead of one 2.5-hour sequence;
`merge_hilbert_derived_8r1.py` reassembles them and re-derives the global gates
from the rows rather than trusting each shard's summary.

### What the derived remainders show

The refusal is not a threshold artefact. Γ's certified remainder degrades
smoothly as a node approaches the upper endpoint, measured on panel 1 with a
deliberately generous radius ladder so the limit cannot be blamed on the
ladder's coarseness:

| 1 − x | Γ remainder |
|---|---|
| 3.0e−14 | 1.9e−16 |
| 7.5e−16 | 2.6e−16 |
| 3.0e−16 | 7.6e−13 |
| 1.0e−16 | 5.6e−11 |
| 3.0e−17 | 5.7e−06 |
| 1.0e−17 | 4.8e−03 |
| 5.8e−18 | 5.5e+00 |
| ≤ 5.2e−18 | construction fails outright |

`|y'|` holds at 1.365 across the whole range, so the division by `y'` is not
the cause — the loss is entirely in Γ. Below 1 − x ≈ 5e−18 the failure changes
character: first `KoenigsError: Koenigs step ratio leaves the principal sector`,
then `polynomial part is not bounded away from zero on the disc`. The
trajectory's grading reaches `tmin = 1e-18` by design, so the grid deliberately
places nodes where the construction has nothing left to say.

**Precision does not fix this.** The state data is stored to 64 significant
digits, so it is not a data limit either. Recomputing the admissibility ratio
at 448, 896 and 1792 bits gives

```
node 850:  Tail/|c_0| = 3.5531e-01   at all three
node 862:  Tail/|c_0| = 2.2029e+02   at all three
```

identical to five digits while the ball radius falls from 1e−132 to exactly 0.
The ratio is a property of the function, not of the computation.

### Why only the upper end

The two punctures are not symmetric, and one line of the state model says why:

```
T(x=0, β) = 1.000000000000000000     at every β on the panel
T(x=1, β) = b(β)                     2.066057 / 2.065290 / 2.064524
```

An earlier revision of this file concluded from those two lines that the upper
puncture "moves with β" and the lower one does not. That reads the symptom
correctly and the cause wrongly. `sexp_β(−1) = 0` holds for **every** β, so
`w₂ = sexp(x−2)` vanishes at `x = 1` for every β and the puncture is
β-stationary at both ends. The asymmetry is not in the mathematics.

It is in the data. The state is interpolated in β from seven Lobatto nodes,
and `η(β) = T(1,β) − b(β)` is that interpolation's error: **1e−64 at the panel
edges** — which *are* Lobatto nodes, where the model reproduces its input to
the input's own 64 digits — and 1.5e−18 to 5.0e−18 at the centres. It is not
a Chebyshev truncation; those coefficients decay to 9.7e−49. The error
displaces the zero of `w₂` to

```
    1 − x* = η / T′(1) = 1.510909e−18 / 1.367028 = 1.105251e−18
```

and quadrature node 862 sits at `1 − x = 1.113721e−18`. The displaced zero
lands on the node to three digits — exactly where the worst ratio was
measured. Below the state model's own β-accuracy the trajectory does not
determine the trace.

### The endpoint model

[`docs/UPPER_ENDPOINT_MODEL.md`](docs/UPPER_ENDPOINT_MODEL.md) derives what
belongs at `x = 1`. In brief: Γ is the Koenigs–Abel velocity `∂_βα/∂_yα`, it
satisfies

```
    Γ(x) = c·T(x)·[ Γ(x−1) − T(x−1) ] ,
```

and `sexp_b(−1) = 0` for every `b` gives `Γ(0) = 0` — the lower endpoint model,
derived rather than asserted — and one step further

```
    Γ(1, β) = −c(β)·b(β) = −e^{β + e^β} ,
```

exact and entire in β, confirmed against the orbit to twelve digits. The
relation needs only **one** logarithm, of an argument near `b`; the step that
kills the orbit is the second one, whose argument goes to zero.

Applied to the 288 rows with `1 − x ≤ 1e−16` (48 per panel, index 816…863,
worst charge 1.50e−13), all 172 refusals disappear and the chain closes. The
cost is bounded by the grading itself: those rows carry a combined quadrature
weight of 9.90e−17 out of 1, and the cotangent kernel against the uniform RH
grid is under 31, so the induced change in the Hilbert integral stays below
4e−28 — against an accepted `HILBERT_REMAINDER = 2e-10`.

That argument is also **tested against the result**: inflating exactly those
288 rows' remainders by 10⁶ and re-running the replay moves `Y` by 2.31e−4
relative (`3.15155307247495e-9 → 3.152281803444456e-9`). The certificate's
dependence on the endpoint model is thus ~2.3e−10 relative. The test could
have refuted the weight argument and did not.

## A defect in the closure's tail rule, and its fix

With both interval producers available, their enclosures can be compared
*as values*: every row, every array, evaluated at five real base samples of
the panel with ball arithmetic, and the two balls must overlap. Coefficient-wise
agreement is not claimed by either producer — both fold their remainder into
the constant coefficient — so that is the right comparison, and it is the one
`verify_primitives_8r1.py` also makes for Γ and the tail.

Eleven of the twelve arrays agree everywhere. `tail_partial_T`, the derivative
of the reciprocal-product tail with respect to the state seed, does **not**:
254 of 2880 sample values have disjoint enclosures, worst at panel 3, row 5,
`s = 1`, where the closure claimed `−0.264981282 ± 8.5e−10` and the constant
producer `−0.264981290520 ± 3.7e−13`. The constant producer sums one more
tower level explicitly; the closure bounds it. Computing that level for the
row in question settles who is right:

```
first neglected d/dT term      8.2885e-9
closure's claimed remainder    6.1562e-10       (factor 3 too small)
value term / its remainder     1.1228e-10 / 1.1423e-10   (fine)
d/dT' term / its remainder     1.2328e-10 / 1.2542e-10   (fine)
```

The cause is in the derivation. `tails.py` argued that "the same ratio
identity applied to the differentiated recursions" makes both derivative
tails geometric with the value ratio `1/(c y_{j+1})`. That is true for the
`T′` seed, on which the tower does not depend at all, and false for the `T`
seed: with `∂_y y_{j+1} = c y_{j+1} ∂_y y_j` the ratio is

```
    ∂_y a_{j+1} / ∂_y a_j  =  (1 + c y_j θ_{j−1}) / (c y_{j+1}),   θ ∈ [0, 1],
```

a factor `1 + c y_j` larger at worst — about twenty at the truncation index of
this segment. The corrected rule bounds the derivative ratio by
`(1 + u_N)/(c e^{u_N})` at the certified lower bound `u_N` of `c y_N`, using
that `(1+u)e^{−u}` decreases and that the tower is monotone above the
threshold `log(2/c)` (which needs `log(2/c) < 2`, now checked next to the two
invariance premises). At the worst row the corrected remainder is `1.2369e−8`
against the neglected `8.2885e−9`.

`test_tails.py` keeps this as a regression guard: at four rows it computes the
first neglected term of the value and of both sensitivities explicitly, one
tower level past the certified index, pointwise at five base samples, and
requires each to lie below the claimed remainder; at the row above it also
asserts that the *old* rule fails that check.

Two remarks on what this means. First, the closure document's Satz 3 is about
the value tail and is unaffected; the derivative claim lived only in the
module. Second, the constant producer is not thereby vindicated: it takes the
next tower level as a β-Taylor series whose spread across the panel exceeds
one, charges `2e-16` for a truncation it never checks, and the closure refuses
to model that level for exactly that reason. Its arrays happened to be right
here because the true term is small; the closure's bound was wrong because it
was derived wrongly. The fix makes the derived `tail_partial_T` remainders
wider by that factor of about twenty on the affected rows and moves nothing
visible in the certificate.

## Certificate values, and a discrepancy with the paper

The chain closes at the declared order. `replay_8r1.py` with both legs
derived (2026-09-06) gives

```
Y    = [ 3.151553145348560e-9 +/- 8.67e-25]
L    = [13613.30588114737     +/- 8.10e-12]
q    = [ 0.1495393411945083   +/- 7.33e-17]
r    = [ 4.63212659033378e-9  +/- 4.40e-24]
p(r) = [-7.8788828633714e-10  +/- 1.86e-24]
Z1 … Z7: PASS
```

with `q < 1` and `Y + (q−1)·r = −7.8788828633714e−10 < 0`. The mixed-provenance
chain (interval leg on the constant, Hilbert leg derived; 2026-08-05, report
`segment_8r1_replay_derived_report.json`) gives `Y = 3.15155307247495e−9`,
`L = 13613.30588130968`, `q = 0.1495393411962912`, `r = 4.63212648323446e−9`,
`p(r) = −7.8788826811874e−10` — the same numbers to eight digits. Routing the
interval leg through the closure tightened the Γ and Γ′ enclosures by two
orders (their widest constant term went from 1.8e−12 to 1.8e−14 and from
3.8e−10 to 1.8e−11) and widened the tails from 1.3e−15 to 1.1e−10 where the
closure refuses to model a further tower level; the certificate barely notices
either, because `Y` is dominated by the Hilbert-side residual and by `Y_tail`.

**Against the five values quoted in Paper VI**, which come from the order-10
computation:

| | Paper VI (order 10) | here (order 22, all derived) | factor |
|---|---|---|---|
| Y | 2.89398623506807e−10 | 3.151553145348560e−9 | **10.890** |
| r | 4.2535568887144e−10 | 4.63212659033378e−9 | **10.890** |
| p(r) | −7.234965587670e−11 | −7.8788828633714e−10 | **10.890** |
| L | 13613.30588130968 | 13613.30588114737 | 1.000 |
| q | 0.1495393411962912 | 0.1495393411945083 | 1.000 |

The three β-dependent quantities scale by *the same* factor, which is what the
structure requires: `r` follows from `Y`, and `p(r) = Y + (q−1)r`, so one
change propagates through all three. `L` and `q` depend on neither the β
truncation nor the Hilbert term and move only in the eleventh digit, where
the tighter Γ′ enclosures enter the Lipschitz majorant.

So the published residual bound is about eleven times tighter than what this
chain can justify. Two things went into that gap and they pull the same way:
the order-10 run declared `"beta_series_order": 22` while coefficients 10…21
were identically zero, so its error budget omitted terms it assumed present;
and its Hilbert rows carried a set constant rather than a derived remainder.
**Paper VI's table in §0 needs the corrected numbers.** The certificate's
*conclusion* is unaffected — the contraction still closes, with `p(r)` negative
by a wide margin — but the stated intervals are not reproducible here.

An earlier revision of this section reported `Y = 2.89398623507837e−10` from a
chain whose Hilbert primitives were still truncated at order 10; that number
is superseded by the one above and is kept out of the table to avoid three
competing values for the same quantity.

## Independent verification

`verify_primitives_8r1.py` regenerates the Koenigs and tail primitives from
the raw trajectory and compares them with the constant producer's arrays. Run
at the declared order 22 on 2026-09-06 (1715 s, concurrently with other work),
verdict in `segment_8r1_verification_report.json`:

```
independent primitive verification: PASS
  worst gamma_remainder          1.6983494e-14
  worst gamma_prime_remainder    1.7815514e-11
  worst tail_remainder           1.1423317e-10
  worst log_kappa_error          2.5006388e-17
  worst base_error               1.3128767e-14
  worst tail_ratio               7.1194233e-9
  phase margin lower             1.8440688
  max |Gamma_new - Gamma_stored|  0
  max |tail_new - tail_stored|    1.1228341e-10
```

How to read the two comparison lines. Both are `.abs_lower()` distances
between the verifier's *polynomial part* at five base samples and the stored
ball, so `0` for Γ means the enclosures are **consistent** — the difference
ball contains zero — not that the arrays are bit-identical; measured with
`.abs_upper()` the two Γ balls sit up to 2.3e−12 apart, and that width is the
constant producer's: its balls carry radii around 1.8e−12 after ~195 recursion
steps, four orders wider than the `BETA_SERIES_REMAINDER = 2e-16` its payload
declares, while the closure models come back at 1e−104 plus a derived
remainder of at most 1.7e−14. The tail line is *not* an inconsistency either:
the verifier's polynomial part stops one tower level earlier than the constant
producer, and its certified remainder (1.14e−10) covers the 1.12e−10 gap. The
value-wise comparison of § *A defect in the closure's tail rule* — which
includes each side's remainder — is the stricter test, and after the fix it
finds no disjoint enclosure in 34560 samples.

What this verification is and is not. It shares no code with the *constant*
producer, and that is the comparison it makes. It shares the closure modules
with `certify_derived_8r1.py`, so it does not independently check the derived
payload; for that leg the independent check is the constant producer, in the
other direction. Two producers with disjoint remainder logic, agreeing as
values, is the evidence this directory offers.

Full chain: produce 17 min, derived interval primitives 5 min on six cores
(about 25 min on one), constant interval primitives 7 min, Hilbert 23 min on
six cores, replay 17 s, verify 29 min — under an hour from nothing to a
verified segment certificate on a six-core machine.

Committed from all of this: the trajectory (1.44 MB), the two replay verdicts
(`segment_8r1_replay_all_derived_report.json` is the certificate,
`segment_8r1_replay_derived_report.json` the mixed-provenance cross-check) and
the verification verdict (`segment_8r1_verification_report.json`). The large
arrays (13–20 MB each) are not; `SHA256SUMS.txt` binds what is here to what
regenerates.

## Layout

| Path | Role |
|---|---|
| `produce_8r1.py` | trajectory from the repo's engine; `save_traj_8r1.py` is the historical entry-point name |
| `certify_derived_8r1.py`, `certify_hilbert_derived_8r1.py`, `merge_hilbert_derived_8r1.py` | the proof path: interval and Hilbert/trace primitives with every remainder derived by the closure layer |
| `certify_8r1.py`, `certify_hilbert_8r1.py` | the constant-remainder producers (`BETA_SERIES_REMAINDER = 2e-16`); the interval one is kept as the cross-check against the derived producer, the Hilbert one does not run at the declared order |
| `replay_8r1.py` | independent re-evaluation of the stored arrays; accepts both interval schemas and both Hilbert schemas, and stops before the arithmetic when a payload carries unproved rows |
| `tmodel.py`, `tmodel2.py`, `states.py`, `koenigs.py`, `tails.py` | the closure layer: Taylor models with certified remainders, replacing set constants by derived bounds |
| `verify_primitives_8r1.py` | regenerates Koenigs and tail primitives from the raw trajectory and compares them with the constant producer's arrays |
| `test_tmodel.py`, `test_tmodel2.py`, `test_tails.py` | soundness tests including negative controls; run them directly, they are scripts |
| `probe_field.py` | diagnostic only, **not** on the proof path |
| `docs/` | the paper's own specification, certification report, review resolution and closure document, plus the endpoint-model note written here |
| `../certificates/` | the schema-v2 certificate template the auditor reads, with its own README |

## What is not closed

The closure document is explicit that two of five points remain open: the
graded Gauss–Legendre bound is derived but its complex tube supremum is not
computed, and the microsegment point turned out to be a substantive gap
rather than a numerical one. The uniform proof subdivision comprises 262144
mathematical microsegments, a consequence of the deliberately coarse
full-ball Lipschitz majorant.

Beyond `8r1` the march segments, the segment-uniform end gluing and the cut
continuation are analytically specified but numerically open. This directory
certifies **one** segment. `fatou_backend.certification` reflects that: it
refuses to report a final certified status, by design, until array-level
replay exists for the whole chain.

**The certificate template** listed by the implementation specification,
`research/certificates/paper_vi_certificate_template.json` (schema v2), shipped
in neither package. It now exists, written *from the auditor* rather than from
any computed result: every number and hash is a `REQUIRED` placeholder and
every claim flag is `false`, so `tetration-cert-audit` rejects it and lists
what a producer has to supply. Its key tree is bound to the auditor's own
certificate fixture by a test. The `8r1` replay report is the material for one
`segments[]` entry; a passing certificate additionally needs the gapless chain
from β = 0 and the array-replay gate, which this release keeps fail-closed.
See `research/certificates/README.md`.

**Open items, by weight**, for whoever continues:

1. **Microsegment induction** (closure §8). The step from the cell radii
   polynomial to the segment endpoint bound `endpoint_total` sums 262144 cells
   linearly with a single Grönwall factor `M1_UPPER = 0.01`, and the closure
   shows that constant does not exist in that form. This is the substantive
   gap; it is the paper's, not the code's.
2. **The graded Gauss–Legendre remainder** `HILBERT_REMAINDER = 2e-10` and the
   tube inflation `COMPLEX_DERIVATIVE_INFLATION = 16` are inherited. The
   closure derives the bound (§6) and has the two-variable Taylor models to
   compute its complex tube supremum, but the per-panel ellipse geometry, the
   tube cover and the summation are not implemented. The closure estimates the
   computation at about an hour and the implementation at considerably more.
   Its sensitivity analysis says even a bound three orders worse would not tip
   the certificate, so this is a matter of rigour, not of margin. Note also
   that `M_Γ′` in the replay is a maximum over the 96 grid rows, and the field
   probe of closure §6.2 shows `Γ` behaves like `t·(A + B log t)` at both ends,
   so its derivative is unbounded there; the ×16 inflation is what covers that
   today.
3. **`|Γ(x)| ≤ 1000·x` near the lower endpoint** is measured, not proved
   (`docs/UPPER_ENDPOINT_MODEL.md` §7). The certificate depends on the endpoint
   model at the 2.3e−10 relative level, so the margin is enormous, but the
   bound still needs a proof — which would have to bound `G` near `y = 1`
   directly.
4. **The β-interpolation floor.** The state is interpolated from seven Lobatto
   nodes in β with an error of 1.5e−18 … 6.2e−18 at panel centres, and the
   graded grid reaches `tmin = 1e-18`. The endpoint model handles the rows
   below that floor, but the trajectory itself does not determine them. Either
   the grading stops at the floor or the state carries more β nodes; both
   change the candidate and therefore the paper's construction, so they are
   design decisions for the paper.
5. **Paper VI §0** still quotes the order-10 intervals; `Y`, `r` and `p(r)`
   need the factor 10.890 (table above).
