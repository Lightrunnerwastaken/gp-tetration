# The upper endpoint of the graded trace

`certify_8r1.py` special-cases the lower endpoint of the trace integral with

```python
if x_index == 0:
    gamma = acb_series([], BETA_SERIES_ORDER)      # Gamma = 0
    koenigs_trace = {"endpoint_stopped_model": True, ...}
```

and asserts, in a comment, that `T(0) = 1` is a removable endpoint of the
Koenigs coordinate. Nothing analogous exists at `x = 1`, and the graded
quadrature — whose whole purpose is to cluster nodes at both ends — puts its
last nodes at `1 - x ≈ 1e-18`, where the recursion has nothing left to say.

This note derives what belongs there. The short answer:

> **Γ(1, β) = −c(β)·b(β)**, exactly, with `c = e^β = log b` and `b = e^c`.

It is the exact counterpart of `Γ(0, β) = 0`, it follows from the same fact
about tetration, and it agrees with the orbit computation to twelve digits.

## 1. What Γ is

Reading the recursion in `koenigs_models` back to its definitions, with
`σ` the Koenigs function of `log_b` at the repelling pair and
`λ = c·L` the multiplier:

```
log_s      = n·log λ + log δ_n            = log σ(y)
chi        = log_s / log λ                = α(y),  the Abel function
numerator  = ∂_β log σ − α·∂_β log λ      = log λ · ∂_β α
slog_deriv = s_n/δ_n = ∂_y log σ          = log λ · ∂_y α
gamma      = numerator / slog_deriv       = ∂_β α / ∂_y α
```

so

```
    Γ = ∂_β α / ∂_y α .
```

Since `dα = ∂_y α · dy + ∂_β α · dβ`, the level sets of `α` move with velocity
`dy/dβ = −Γ`. **Γ is minus the velocity of a point under the base flow that
holds the Koenigs–Abel coordinate fixed.**

A caution that matters for reading the rest of the certificate: `α` is the
*Koenigs* Abel function, not Kneser's. They differ by exactly the θ-correction
that the Riemann–Hilbert problem determines, so `Γ ≠ −∂_β T` in general —
measured, the two disagree by 10 % at mid-panel and agree only as `x → 1`.
Γ is the raw velocity the RH step corrects, not the corrected one.

## 2. The functional equation

`α` satisfies `α(log_b y) = α(y) − 1`. Write `w = log_b y`. Differentiating in
β at fixed y, and using `∂_β(log y / c) = −log y / c = −w`:

```
    ∂_β α(w) − w·∂_y α(w) = ∂_β α(y)                                  (A)
```

Differentiating in y, with `∂_y w = 1/(c y)`:

```
    ∂_y α(w) = c·y·∂_y α(y)                                           (B)
```

Divide (A) by `∂_y α(y)` and substitute (B):

```
    G(y) = c·y·[ G(log_b y) − log_b y ]                               (★)
```

where `G(y) := ∂_β α(y)/∂_y α(y)`. With `Γ(x) := G(T(x))` and
`T(x−1) = log_b T(x)` this is

```
    Γ(x) = c·T(x)·[ Γ(x−1) − T(x−1) ] .                               (★′)
```

## 3. Both endpoints, from one fact

**`sexp_b(−1) = 0` for every base b.** The point `y = 0` therefore does not
move under the base flow, so `G(0) = 0`.

Now read (★) at the two relevant points. At `y = 1`, where `log_b 1 = 0`:

```
    G(1) = c·1·[ G(0) − 0 ] = c·G(0) = 0 ,
```

which is `Γ(0) = G(T(0)) = G(1) = 0` — the lower endpoint model, derived
rather than asserted. At `y = b`, where `log_b b = 1`:

```
    G(b) = c·b·[ G(1) − 1 ] = −c·b ,
```

that is

```
    Γ(1, β) = −c(β)·b(β) = −e^{β + e^β} ,
```

an entire function of β. The two endpoint models are the same statement
evaluated one step apart.

**Check.** At `1 − x = 5.2e−15` the orbit gives `Γ = −1.4978944341689227`;
`−c·b` at the same β is `−1.497894434169587674`. Twelve digits.

## 4. Why the lower endpoint computes and the upper one does not

Both ends are singular for the recursion, one step apart:

| | state | first zero | dies at |
|---|---|---|---|
| x → 0 | y → 1 | `w₁ = sexp(x−1) → 0` | step 2 |
| x → 1 | y → b | `w₂ = sexp(x−2) → 0` | step 3 |

Mathematically both zeros sit at a β-**independent** place — `sexp_β(−1) = 0`
for all β — so `log w = log(x − x₀) + log A(β)` separates in both cases and
neither should trouble a β-Taylor model. The asymmetry is not in the
mathematics. It is in the data.

The state is interpolated in β from **seven Lobatto nodes**, and that
interpolation has an error. Measured as `η(β) = T(1,β) − b(β)`:

| panel | t = −1 | t = 0 | t = +1 |
|---|---|---|---|
| 1 | 5.5e−65 | 1.51e−18 | 3.9e−64 |
| 3 | 3.0e−64 | 5.04e−18 | 2.5e−64 |
| 6 | 4.3e−64 | 1.50e−18 | 2.9e−65 |

`η` vanishes at the panel edges — those *are* Lobatto nodes, and there the
model reproduces its input to the input's own 64 digits — and peaks at the
centres. It is the degree-6 interpolation error, not a truncation of the
48-term Chebyshev state (whose coefficients decay to 9.7e−49).

That error displaces the zero of `w₂` from `x = 1` to

```
    1 − x* = η / T′(1) = 1.510909e−18 / 1.367028 = 1.105251e−18 ,
```

and quadrature node 862 sits at `1 − x = 1.113721e−18`. **The displaced zero
lands on the node, to three digits** — which is exactly where the worst
failure was measured (`Σ|c_k|/|c₀| = 220`, against ~1 at the neighbours).

Since `η` depends on β, `x*(β)` sweeps across the panel, and the β-Taylor
model of `w₂` there really does straddle a zero. Below the state model's own
β-accuracy the trajectory does not determine the trace. No amount of orbit
work, precision or radius-ladder refinement recovers those rows: recomputing
the admissibility ratio at 448, 896 and 1792 bits leaves it identical to five
digits while the ball radius falls to exactly 0.

## 5. The model

Evaluate (★′) directly. The point is that it needs only **one** logarithm,
`T(x−1) = log_b T(x)`, whose argument sits near `b` — the step that destroys
the orbit is the *second* one, whose argument goes to zero.

```python
w      = y.log() / c                    # T(x-1), well conditioned
factor = c * y
gamma  = factor * (-w)                  # + factor * Gamma(x-1)
charge = factor.bound() * slope * eps   # |Gamma(x-1)| <= slope * eps
```

`Γ(x−1)` is enclosed by the lower endpoint bound. Measured on the graded grid,
`|Γ(x)|/x` is 861 at `x = 1.02e−18` and falls to 453 at `x = 8e−12` — the same
`≈ 29·log(1/x)` law that governs the upper endpoint slope (353 → 1007 over the
mirrored range). `slope = 1000` covers the whole relevant range with margin.

At `x = 1` the model returns `c·b·(0 − 1) = −c·b`, the exact value of §3.

**Cross-check where both methods work** (panel 1, `ε = 1e−16 … 1e−11`):

| 1 − x | orbit | endpoint model | charge | separation |
|---|---|---|---|---|
| 9.19e−12 | −1.4978944335123443 | −1.4978944341479132 | 1.38e−08 | 6.04e−09 |
| 1.51e−13 | −1.4978944341543836 | −1.4978944341692313 | 2.27e−10 | 1.22e−10 |
| 5.21e−15 | −1.4978944341689227 | −1.4978944341695755 | 7.81e−12 | 4.86e−12 |
| 1.02e−16 | −1.4978944341695712 | −1.4978944341695875 | 1.53e−13 | 1.11e−13 |

The enclosures are consistent at every node, and the separation stays under
the charge throughout. The model is *wider* than the orbit above
`1 − x ≈ 1e−16`, which is why the cut-over is placed there and not higher.

## 6. What it costs

The graded weights of the affected rows are tiny — that is what grading does.
The model applies where `1 − x ≤ 1e-16`, which is 48 of the 864 nodes
(index 816…863), and

```
    Σ |w_t| over the endpoint set  =  9.900e−17      (total weight 1)
```

The RH evaluation points are `k/96`, so the closest one to `t ≈ 1` is `95/96`
and the cotangent kernel is bounded by `cot(π/96) < 31`. With a worst endpoint
charge of `1.496e−13` and `|T′| ≈ 1.37`, the induced change in the Hilbert
integral is bounded by

```
    9.900e−17 · (1.496e−13 / 1.37) · 31  <  4e−28 ,
```

against an accepted `HILBERT_REMAINDER = 2e-10`.

**Measured end to end, not only bounded.** An argument about the size of a
contribution should be tested against the result, not only against its own
derivation. Inflating the remainder of exactly those 288 rows (48 per panel)
by a factor of 10⁶ and re-running the replay moves the residual bound from

```
    Y = 3.15155307247495e-9   →   Y = 3.152281803444456e-9
```

a relative change of 2.31e−4 for a millionfold widening. The certificate's
dependence on the endpoint model is therefore about 2.3e−10 in relative terms,
or 7.3e−19 absolute in Y. The test could have refuted the weight argument; it
did not, and it moves by a proportionate amount rather than not at all, which
is what a live dependency of the right size looks like.

An earlier revision of this section quoted `4.258e−17` — that is the weight sum
over nodes 830…863, the *failing* set, not the larger set the model is applied
to. The corrected figure is 2.3× larger and changes nothing downstream.

## 7. What is proved and what is not

**Proved.** The functional equation (★) and (★′); the two endpoint values
`Γ(0) = 0` and `Γ(1) = −c·b` given `G(0) = 0`; the conditioning claim that
(★′) needs only a logarithm of an argument near `b`.

**Measured, not proved.** The bound `|Γ(x)| ≤ 1000·x` near the lower endpoint
is an interval evaluation on the graded grid, extended to the continuum by the
observed `29·log(1/x)` law. A proof would bound `G` near `y = 1` directly;
the certificate currently assumes something *stronger* without saying so —
`certify_8r1.py` sets `Γ = 0` at `x_index == 0` with no radius at all, i.e. it
takes `slope = 0`.

**Assumed by inheritance.** `HILBERT_REMAINDER = 2e-10` is carried over
unchanged; the closure lists its graded Gauss–Legendre bound as derived with
the complex tube supremum not computed.

**Not addressed here.** The β-interpolation error `η ≈ 1.5e−18 … 6.2e−18` is
the state model's accuracy floor, and it is now the binding constraint on how
finely the graded grid can be resolved. Placing quadrature nodes below it —
as `tmin = 1e-18` does — asks the trajectory for information it does not
carry. Either the grading should stop at `η`, or the state should be
interpolated from more β nodes.
