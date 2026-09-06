"""Base- and state-model construction shared by the closure programs.

The trajectory stores, for each of the seven Chebyshev--Lobatto base nodes,
the Chebyshev coefficients of the state ``T(., x)`` and of ``T'`` on
``x in [0,1]``.  On one base panel the base dependence is carried as a
Taylor model in the scaled panel variable ``s``, where

    beta = beta_mid + beta_half * s ,     |s| <= 1 .

Because the seven stored nodes determine a degree-six interpolant, the base
models are *polynomials*: they carry no truncation remainder at all, and the
Taylor-model remainder of every derived quantity therefore originates only in
the transcendental steps.
"""

from __future__ import annotations

from flint import acb, arb

from tmodel import TM


def cheb_coefficients(values: list) -> list:
    """Chebyshev--Lobatto coefficients of the interpolant through ``values``."""

    n = len(values) - 1
    coefficients = []
    for k in range(n + 1):
        total = values[0] / 2 + values[-1] * ((-1) ** k) / 2
        for j in range(1, n):
            total += values[j] * (arb.pi() * j * k / n).cos()
        coefficient = 2 * total / n
        if k in (0, n):
            coefficient /= 2
        coefficients.append(coefficient)
    return coefficients


def cheb_eval_tm(coefficients: list[TM], point: TM, order: int) -> TM:
    """Clenshaw evaluation of a Chebyshev series of Taylor models."""

    b1 = TM.zero(order)
    b2 = TM.zero(order)
    for coefficient in coefficients[:0:-1]:
        b0 = point * b1 * 2 - b2 + coefficient
        b2, b1 = b1, b0
    return point * b1 - b2 + coefficients[0]


def cheb_eval_acb(coefficients: list[acb], point: acb) -> acb:
    b1 = acb(0)
    b2 = acb(0)
    for coefficient in coefficients[:0:-1]:
        b0 = 2 * point * b1 - b2 + coefficient
        b2, b1 = b1, b0
    return point * b1 - b2 + coefficients[0]


def state_second_derivative(derivative_coefficients: list[arb]) -> list[arb]:
    """Chebyshev coefficients of ``T''`` from those of ``T'``."""

    n = len(derivative_coefficients) - 1
    output = [arb(0) for _ in range(n + 1)]
    if n == 0:
        return output
    output[n - 1] = 2 * n * derivative_coefficients[n]
    for k in range(n - 2, -1, -1):
        output[k] = 2 * (k + 1) * derivative_coefficients[k + 1]
        if k + 2 <= n:
            output[k] += output[k + 2]
    output[0] /= 2
    return [2 * value for value in output]


class PanelStates:
    """Chebyshev state models on one base panel, as Taylor models in ``s``."""

    def __init__(self, trajectory: dict, panel_index: int, order: int):
        self.order = order
        self.panel_index = panel_index
        beta_a = arb(trajectory["segment"]["beta_a"])
        beta_b = arb(trajectory["segment"]["beta_b"])
        nodes = [
            arb(value)
            for value in trajectory["transforms"]["base_lobatto_beta_nodes"]
        ]
        self.beta_left = nodes[panel_index]
        self.beta_right = nodes[panel_index + 1]
        self.beta_mid = (self.beta_left + self.beta_right) / 2
        self.beta_half = (self.beta_right - self.beta_left) / 2
        self.beta = TM.linear(
            acb(self.beta_mid), acb(self.beta_half), order
        )
        base_point = (self.beta - beta_a) * (2 / (beta_b - beta_a)) - 1

        coefficients = [
            [arb(value) for value in node["state"]["coefficients"]]
            for node in trajectory["nodes"]
        ]
        derivatives = [
            [arb(value) for value in node["state"]["derivative_coefficients"]]
            for node in trajectory["nodes"]
        ]
        seconds = [state_second_derivative(row) for row in derivatives]

        self.state = self._models(coefficients, base_point)
        self.state_prime = self._models(derivatives, base_point)
        self.state_second = self._models(seconds, base_point)

    def _models(self, per_node: list[list[arb]], base_point: TM) -> list[TM]:
        count = len(per_node[0])
        models = []
        for index in range(count):
            # Stored base nodes ascend; the Lobatto transform is descending.
            values = [per_node[node][index] for node in range(7)][::-1]
            coefficients = [
                TM.constant(acb(value), self.order)
                for value in cheb_coefficients(values)
            ]
            models.append(
                cheb_eval_tm(coefficients, base_point, self.order)
            )
        return models

    def at(self, x) -> tuple[TM, TM, TM]:
        """State, first and second ``x``-derivative at a (possibly complex) x."""

        point = TM.constant(acb(x) * 2 - 1, self.order)
        return (
            cheb_eval_tm(self.state, point, self.order),
            cheb_eval_tm(self.state_prime, point, self.order),
            cheb_eval_tm(self.state_second, point, self.order),
        )
