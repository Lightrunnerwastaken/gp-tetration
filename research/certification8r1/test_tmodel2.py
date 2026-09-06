"""Soundness tests for the two-variable Taylor models."""

from __future__ import annotations

from flint import acb, arb, ctx

from tmodel import TM, configure
from tmodel2 import TM2


NS = 10
NT = 8


def build_samples() -> list[tuple[acb, acb]]:
    root = (arb(1) / 2).sqrt()
    points = [
        acb(1),
        acb(-1),
        acb("0.5"),
        acb(0, 1),
        acb(root, root),
        acb(0),
    ]
    return [(s, t) for s in points for t in points]


def check(name: str, model: TM2, truth, s: acb, t: acb) -> bool:
    value = acb(0)
    power = acb(1)
    for coefficient in model.a:
        value += coefficient.poly()(s) * power
        power = power * t
    # Containment test: the true value must lie in the evaluated coefficient
    # ball inflated by the remainder.  ``abs_lower`` of the difference is the
    # smallest modulus attained on that ball, so the test is exactly
    # ``truth in value + D(0, r)``.
    difference = (truth - value).abs_lower()
    ok = bool(difference <= model.r)
    if not ok:
        print(
            f"  [FAIL] {name}: |err|={difference.str(6)} r={model.r.str(6)}"
        )
    return ok


def main() -> int:
    configure(NS, 448)
    samples = build_samples()
    failures = 0

    beta_mid = arb("-0.3283344913544076")
    beta_half = arb("1.9089214613628350e-3")
    t_mid = arb("0.3")
    t_half = arb("0.03")

    beta = TM.linear(acb(beta_mid), acb(beta_half), NS)
    tau = TM2.linear(TM.constant(acb(t_mid), NS), TM.constant(acb(t_half), NS), NT)
    base = TM2.constant(beta, NT)

    combined = base + tau * tau + tau * base
    exponential = combined.exp()
    logarithm = (combined + 4).log()
    inverse = (combined + 4).inverse()
    quotient = exponential / (combined + 4)

    for s, t in samples:
        b = acb(beta_mid) + acb(beta_half) * s
        x = acb(t_mid) + acb(t_half) * t
        raw = b + x * x + x * b
        failures += not check("combined", combined, raw, s, t)
        failures += not check("exp", exponential, raw.exp(), s, t)
        failures += not check("log", logarithm, (raw + 4).log(), s, t)
        failures += not check("inverse", inverse, acb(1) / (raw + 4), s, t)
        failures += not check(
            "quotient", quotient, raw.exp() / (raw + 4), s, t
        )

    print(f"orders: n_s = {NS}, n_tau = {NT}, precision = {ctx.prec} bits")
    for name, model in (
        ("combined", combined),
        ("exp", exponential),
        ("log", logarithm),
        ("inverse", inverse),
        ("quotient", quotient),
    ):
        print(
            f"  {name:10s} bound = {model.bound().str(8):>16s} "
            f"remainder = {model.r.str(8)}"
        )

    print()
    if failures:
        print(f"{failures} failing checks")
        return 1
    print(f"all TM2 invariants hold on {len(samples)} bidisc samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
