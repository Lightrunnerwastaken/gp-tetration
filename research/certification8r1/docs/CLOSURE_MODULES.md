# `8r1` closure — derived residual budgets

This directory replaces the set constants of the `8r1` pre-certificate by
proved interval bounds, and supplies a primitive verifier that shares no code
with the original producers.

## Files

| File | Role |
|---|---|
| `tmodel.py` | one-variable Taylor models with a certified remainder disc: exact convolution defect for multiplication and inversion, certified Cauchy tails for `exp`/`log`, branch-safe logarithm, a-posteriori fixed-point certificate for `L = exp(cL)` |
| `tmodel2.py` | two-variable Taylor models (a model in `tau` whose coefficients are one-variable models in `s`); needed because a plain interval in the base parameter is destroyed by wrapping within twenty Koenigs steps |
| `states.py` | Chebyshev state models on one base panel; these are polynomials, so they carry no truncation remainder |
| `koenigs.py` | stopped Koenigs coordinate with the derived remainder that replaces the factors `128` and `1024` |
| `tails.py` | reciprocal-product tail with the proved ratio rule that replaces `R <= 2 abs(a_N)`; the `d/dT` sensitivity carries its own, larger ratio `(1 + c y_j theta)/(c y_{j+1})` -- corrected 2026-09-06, see the module docstring and `test_tails.py` |
| `verify_primitives_8r1.py` | independent regeneration of the Koenigs and tail primitives from the raw trajectory, plus comparison against the stored pre-certificate arrays |
| `probe_field.py` | diagnostic probe of the geometric field off the real axis; **not** part of the proof path |
| `test_tmodel.py`, `test_tmodel2.py` | soundness tests: model invariants against independent point evaluation, plus negative controls |
| `test_tails.py` | soundness test for the tail remainders: the first neglected tower level, computed explicitly and pointwise, must lie below each claimed remainder (value, `d/dT`, `d/dT'`); negative control against the pre-2026-09-06 derivative rule |

## Requirements

```bash
pip install --break-system-packages python-flint==0.9.0 mpmath
```

PARI/GP is **not** required: nothing here starts the target-base engine or
reads target-base values. The verifier consumes only `trajectory_8r1.json`.

## Running

```bash
python3 test_tmodel.py
python3 test_tmodel2.py
python3 test_tails.py

python3 verify_primitives_8r1.py \
  --trajectory        path/to/trajectory_8r1.json \
  --stored-primitives path/to/segment_8r1_interval_primitives.json \
  --output            segment_8r1_derived_primitives.json \
  --report            segment_8r1_independent_verification.json
```

## One trap worth naming

`python-flint` truncates every series operation at the global `ctx.cap`,
whose default is `10`. A model built with twenty-two coefficients is
silently cut back unless `ctx.cap` is raised, and the discarded coefficients
are invisible in the serialized output — the pre-certificate declares order
22 and computed order 10. Every entry point here calls
`tmodel.configure(order, bits)` first, and `tmodel.assert_full_length`
re-checks the achieved length after each series operation.
