"""Soundness checks for ``tails.certified_tail`` against one more tower level.

The certified tail stops at an index ``N`` and charges the levels ``j > N`` to
a geometric remainder.  The first neglected term dominates that remainder, so
computing it *explicitly* -- one more level of the very same recursion -- is a
necessary condition the claimed remainders must pass, for the value and for
both seed sensitivities.  The extra level is evaluated pointwise, at real
base samples ``s`` of the panel, with scalar ball arithmetic: as a Taylor
model in the base it is not representable (that is why the closure stops
there), as a number at a fixed base it is harmless.  Run this file directly;
it is a script, like ``test_tmodel.py``.

It also keeps a negative control.  Until 2026-09-06 the module charged the
``d/dy`` sensitivity tail with the *value* ratio ``1/(c y_{N+1})``.  The true
ratio carries the factor ``1 + c y_j theta`` (see the module docstring), and
at the truncation index of this segment that is a factor of about twenty.
The defect was found by comparing the closure's ``tail_partial_T`` enclosures
with those of the constant producer, which sums one more level explicitly:
254 of 2880 sample values did not overlap.  The control below asserts that
the old rule *fails* the explicit check, so the comparison cannot silently
go back to it.
"""

from __future__ import annotations

import json
from pathlib import Path

from flint import acb, arb

from koenigs import KoenigsBase
from states import PanelStates
from tails import certified_tail
from tmodel import TM, configure


ORDER = 22
PRECISION_BITS = 448
TAIL_STOP = "1e-50"
SAMPLES = ("-1", "-0.5", "0", "0.5", "1")
# (panel, x_index): the worst row of the 2026-09-06 comparison, one interior
# row, one row near the lower endpoint and one near the upper endpoint.
ROWS = ((3, 5), (0, 48), (2, 1), (5, 90))


def explicit_terms(c: acb, y: acb, y_prime: acb, count: int):
    """Scalar ``a_j``, ``d_y a_j``, ``d_y' a_j`` for ``j = 1 .. count``.

    The recursion of ``certified_tail`` at one fixed base, without any
    truncation logic, so the caller can look one level past the certified
    stopping index.
    """

    current_y, current_y_prime = y, y_prime
    dy_dy, dy_dyp, dyp_dy, dyp_dyp = acb(1), acb(0), acb(0), acb(1)
    terms, sens_y, sens_yp = [], [], []
    for _ in range(count):
        inverse_prime = 1 / current_y_prime
        terms.append(inverse_prime / c)
        sens_y.append(-(dyp_dy * inverse_prime * inverse_prime / c))
        sens_yp.append(-(dyp_dyp * inverse_prime * inverse_prime / c))
        next_y = (c * current_y).exp()
        next_dy_dy = c * next_y * dy_dy
        next_dy_dyp = c * next_y * dy_dyp
        next_y_prime = c * next_y * current_y_prime
        next_dyp_dy = c * (next_dy_dy * current_y_prime + next_y * dyp_dy)
        next_dyp_dyp = c * (next_dy_dyp * current_y_prime + next_y * dyp_dyp)
        current_y, current_y_prime = next_y, next_y_prime
        dy_dy, dy_dyp = next_dy_dy, next_dy_dyp
        dyp_dy, dyp_dyp = next_dyp_dy, next_dyp_dyp
    return terms, sens_y, sens_yp


def at(model: TM, s: acb) -> acb:
    return model.poly()(s)


def main() -> int:
    configure(ORDER, PRECISION_BITS)
    stop = arb(TAIL_STOP)
    trajectory = json.loads(
        Path(__file__).with_name("trajectory_8r1.json").read_text(
            encoding="utf-8"
        )
    )
    samples = [acb(arb(text)) for text in SAMPLES]
    failures = 0
    control_seen = False
    for panel_index, x_index in ROWS:
        states = PanelStates(trajectory, panel_index, ORDER)
        base = KoenigsBase(states.beta, ORDER)
        y, y_prime, _ = states.at(arb(x_index) / 96)
        tail, tail_y, tail_yp, trace = certified_tail(
            base.c, y, y_prime, stop
        )
        n = trace["terms"] - 1
        worst = {"value": arb(0), "partial_T": arb(0), "partial_T_prime": arb(0)}
        old_worst = arb(0)
        for s in samples:
            terms, sens_y, sens_yp = explicit_terms(
                at(base.c, s), at(y, s), at(y_prime, s), n + 1
            )
            for name, sequence in (
                ("value", terms),
                ("partial_T", sens_y),
                ("partial_T_prime", sens_yp),
            ):
                worst[name] = max(worst[name], sequence[n].abs_upper())
            # what the pre-2026-09-06 rule would have claimed at this sample
            old_worst = max(
                old_worst, sens_y[n - 1].abs_upper() * trace["geometric_factor"]
            )
        claimed = {
            "value": trace["remainder"],
            "partial_T": trace["remainder_partial_T"],
            "partial_T_prime": trace["remainder_partial_T_prime"],
        }
        print(f"panel {panel_index} x_index {x_index}: N = {n}")
        for name in worst:
            ok = bool(worst[name] <= claimed[name])
            failures += not ok
            print(
                f"  [{'ok ' if ok else 'BAD'}] {name:16s} first neglected "
                f"{worst[name].str(5)}  <=  claimed {claimed[name].str(5)}"
            )

        if (panel_index, x_index) == (3, 5):
            # Negative control: the value-ratio rule must NOT cover the
            # neglected d/dy term at the row where the defect was found.
            old_rule_fails = not bool(worst["partial_T"] <= old_worst)
            control_seen = True
            failures += not old_rule_fails
            print(
                f"  [{'ok ' if old_rule_fails else 'BAD'}] negative control: "
                f"value-ratio rule would claim {old_worst.str(5)} against a "
                f"neglected term of {worst['partial_T'].str(5)}"
            )
    if not control_seen:
        failures += 1
        print("negative control row missing")
    print()
    print("tails: " + ("all checks hold" if not failures else f"{failures} FAILED"))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
