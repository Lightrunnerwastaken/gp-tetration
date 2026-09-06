"""Soundness and sharpness tests for the certified Taylor models.

For every model we check the defining invariant at sample points of the
closed unit disc: the true value, computed independently with plain Arb ball
arithmetic at a thin argument, must lie inside ``P(s) + D(0,r)``.
"""

from __future__ import annotations

import sys

from flint import acb, acb_series, arb, ctx

from tmodel import TM, TaylorModelError, configure, fixed_point_exp


def true_fixed_point(c_point: acb, seed: acb) -> acb:
    """Newton refinement of ``L = exp(c L)`` on the branch through ``seed``.

    Using the defining equation rather than a Lambert-``W`` call keeps the
    reference value on the same analytic branch as the model for complex
    samples, where the principal ``W_{-1}`` labelling is discontinuous.
    """

    z = acb(seed)
    for _ in range(80):
        value = (c_point * z).exp()
        z = z - (z - value) / (1 - c_point * value)
    return z


ORDER = 22


def build_samples() -> list[acb]:
    """Sample points of the closed unit disc, built at working precision.

    The list is constructed only after :func:`configure`, because an inexact
    decimal such as ``0.6`` parsed at the default 53-bit precision would be a
    ball of radius ``1e-17``; the comparison would then measure the width of
    the sample rather than the width of the model.
    """

    root = (arb(1) / 2).sqrt()
    return [
        acb(1),
        acb(-1),
        acb("0.5"),
        acb("-0.75"),
        acb(0, 1),
        acb(0, -1),
        acb(root, root),
        acb(-root, root),
        acb(0),
    ]


def check(name: str, model: TM, truth, sample: acb) -> bool:
    value = model.poly()(sample)
    difference = (truth - value).abs_upper()
    ok = bool(difference <= model.r)
    status = "ok " if ok else "FAIL"
    print(
        f"  [{status}] {name:28s} |err|={difference.str(6):>14s} "
        f"r={model.r.str(6):>14s}"
    )
    return ok


def main() -> int:
    configure(ORDER, 448)
    samples = build_samples()
    failures = 0
    print(f"series cap = {ctx.cap}, precision = {ctx.prec} bits")
    print()

    beta_mid = arb("-0.3283344913544076")
    beta_half = arb("1.9089214613628350e-3")
    beta = TM.linear(acb(beta_mid), acb(beta_half), ORDER)

    print("exp / log / inverse invariants")
    c = beta.exp()
    log_c = c.log()
    inv_c = c.inverse()
    quotient = beta / c
    for sample in samples:
        argument = acb(beta_mid) + acb(beta_half) * sample
        failures += not check("exp(beta)", c, argument.exp(), sample)
        failures += not check("log(exp(beta))", log_c, argument, sample)
        failures += not check(
            "1/exp(beta)", inv_c, acb(1) / argument.exp(), sample
        )
        failures += not check(
            "beta/exp(beta)", quotient, argument / argument.exp(), sample
        )

    print()
    print("remainder magnitudes (must be tiny, not a set constant)")
    for name, model in (
        ("exp(beta)", c),
        ("log(exp(beta))", log_c),
        ("1/exp(beta)", inv_c),
        ("beta/exp(beta)", quotient),
    ):
        print(f"  {name:20s} r = {model.r.str(8)}")

    print()
    print("repelling fixed point L = exp(c L)")
    c_series = acb_series(
        [acb(beta_mid), acb(beta_half)], ORDER
    ).exp()
    seed = -(-c_series).lambertw(-1) / c_series
    fixed, kappa = fixed_point_exp(c, seed, ORDER)
    print(f"  contraction factor kappa = {kappa.str(8)}")
    print(f"  certified remainder      = {fixed.r.str(8)}")
    for sample in samples:
        argument = acb(beta_mid) + acb(beta_half) * sample
        c_point = argument.exp()
        truth = true_fixed_point(c_point, fixed.poly()(sample))
        failures += not check("L(beta)", fixed, truth, sample)

    print()
    print("fixed-point equation residual on the whole disc")
    residual = (fixed - (c * fixed).exp()).bound()
    print(f"  sup |L - exp(cL)| <= {residual.str(8)}")

    print()
    print("negative controls (must raise)")
    for name, thunk in (
        (
            "inverse of a model straddling zero",
            lambda: TM([acb(0), acb(1)], arb(0), ORDER).inverse(),
        ),
        (
            "log of a model straddling zero",
            lambda: TM([acb("0.1"), acb(1)], arb(0), ORDER).log(),
        ),
        (
            "inverse whose remainder reaches zero",
            lambda: TM([acb(1)], arb(2), ORDER).inverse(),
        ),
    ):
        try:
            thunk()
        except TaylorModelError as error:
            print(f"  [ok ] {name}: {error}")
        else:
            print(f"  [FAIL] {name} did not raise")
            failures += 1

    print()
    if failures:
        print(f"{failures} failing checks")
        return 1
    print("all Taylor-model invariants hold")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
