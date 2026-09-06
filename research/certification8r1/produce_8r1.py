"""Recompute the missing Paper-VI edge trajectory ``8r1``.

This program is a *candidate producer*.  It deliberately keeps the target
Kneser engine on the producer side.  The independent checker never imports
this module and never starts PARI/GP.

The Riemann--Hilbert field is not differentiated from the target trajectory.
It is rebuilt from the single-base dynamics using the Koenigs recursions of
Paper V and a graded composite Gauss--Legendre Hilbert transform.  Target
engine beta stencils are exported only as an independent numerical diagnostic.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
from pathlib import Path
from typing import Iterable

import mpmath as mp

from fatou_backend.gp_backend import FatouGP


DPS = 80
STATE_N = 48
BASE_N = 7
RH_N = 96
RH_HEAD = 32
PANEL_ORDER = 24
TMIN = mp.mpf("1e-18")
KOENIGS_STOP = mp.mpf("1e-45")
RAW_DEPTH = 4
STENCIL_STEP = mp.mpf("1e-10")


def decimal(value: mp.mpf | mp.mpc, digits: int = 64) -> str:
    if isinstance(value, mp.mpc):
        raise TypeError("use complex_record for complex values")
    return mp.nstr(value, digits, min_fixed=0, max_fixed=0)


def complex_record(value: mp.mpc, digits: int = 64) -> dict[str, str]:
    return {
        "re": decimal(mp.re(value), digits),
        "im": decimal(mp.im(value), digits),
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def cheb_lobatto_nodes(count: int) -> list[mp.mpf]:
    n = count - 1
    return [(mp.cos(mp.pi * j / n) + 1) / 2 for j in range(count)]


def cheb_coefficients(values: list[mp.mpf | mp.mpc]) -> list[mp.mpf | mp.mpc]:
    n = len(values) - 1
    coeffs: list[mp.mpf | mp.mpc] = []
    for k in range(n + 1):
        total = values[0] / 2 + values[-1] * ((-1) ** k) / 2
        for j in range(1, n):
            total += values[j] * mp.cos(mp.pi * j * k / n)
        coeffs.append(2 * total / n)
    coeffs[0] /= 2
    coeffs[-1] /= 2
    return coeffs


def cheb_derivative(coeffs: list[mp.mpf | mp.mpc]) -> list[mp.mpf | mp.mpc]:
    n = len(coeffs) - 1
    deriv: list[mp.mpf | mp.mpc] = [mp.mpf(0)] * (n + 1)
    if n == 0:
        return deriv
    deriv[n - 1] = 2 * n * coeffs[n]
    for k in range(n - 2, -1, -1):
        deriv[k] = 2 * (k + 1) * coeffs[k + 1]
        if k + 2 <= n:
            deriv[k] += deriv[k + 2]
    deriv[0] /= 2
    # The stored polynomial uses s = 2*x-1.
    return [2 * value for value in deriv]


def cheb_eval(coeffs: list[mp.mpf | mp.mpc], x: mp.mpf | mp.mpc):
    s = 2 * x - 1
    b1 = mp.mpc(0) if any(isinstance(c, mp.mpc) for c in coeffs) else mp.mpf(0)
    b2 = b1
    for k in range(len(coeffs) - 1, 0, -1):
        b0 = 2 * s * b1 - b2 + coeffs[k]
        b2, b1 = b1, b0
    return s * b1 - b2 + coeffs[0]


def fourier_coefficients(values: list[mp.mpf | mp.mpc]) -> dict[int, mp.mpc]:
    count = len(values)
    result: dict[int, mp.mpc] = {}
    for k in range(-count // 2 + 1, count // 2):
        total = mp.mpc(0)
        for j, value in enumerate(values):
            total += value * mp.exp(-2j * mp.pi * k * j / count)
        result[k] = total / count
    return result


def barycentric_weights(nodes: list[mp.mpf]) -> list[mp.mpf]:
    n = len(nodes) - 1
    return [
        mp.mpf((-1) ** j) * (mp.mpf("0.5") if j in (0, n) else 1)
        for j in range(n + 1)
    ]


def differentiation_matrix(nodes: list[mp.mpf]) -> list[list[mp.mpf]]:
    weights = barycentric_weights(nodes)
    size = len(nodes)
    matrix = [[mp.mpf(0) for _ in range(size)] for _ in range(size)]
    for i in range(size):
        for j in range(size):
            if i != j:
                matrix[i][j] = weights[j] / (weights[i] * (nodes[i] - nodes[j]))
        matrix[i][i] = -sum(matrix[i][j] for j in range(size) if j != i)
    return matrix


def panel_rule() -> tuple[list[mp.mpf], list[mp.mpf], list[tuple[mp.mpf, mp.mpf]]]:
    rule_nodes, rule_weights = mp.gauss_quadrature(PANEL_ORDER, "legendre")
    left = [TMIN]
    point = TMIN
    while point < mp.mpf("0.1"):
        point = min(point * 10, mp.mpf("0.1"))
        left.append(point)
    left.append(mp.mpf("0.5"))
    breaks = left + [1 - value for value in reversed(left[:-1])]
    panels = list(zip(breaks[:-1], breaks[1:]))
    nodes: list[mp.mpf] = []
    weights: list[mp.mpf] = []
    for a, b in panels:
        midpoint = (a + b) / 2
        halfwidth = (b - a) / 2
        for node, weight in zip(rule_nodes, rule_weights):
            nodes.append(midpoint + halfwidth * node)
            weights.append(halfwidth * weight)
    return nodes, weights, panels


class KoenigsData:
    def __init__(self, beta: mp.mpf):
        self.beta = beta
        self.c = mp.exp(beta)
        self.base = mp.exp(self.c)
        self.L = -mp.lambertw(-self.c, -1) / self.c
        self.lam = self.c * self.L
        self.log_lam = mp.log(self.lam)
        self.dL = self.c * self.L * self.L / (1 - self.lam)
        self.dlam = self.c * (self.dL + self.L)
        self.dlog_lam = self.dlam / self.lam

    def geometric_field(self, y: mp.mpf, y_prime: mp.mpf) -> tuple[mp.mpc, dict]:
        """Return A'(y) Gamma(y), using the exact Paper-V recursions."""

        w = mp.mpc(y)
        v = mp.mpc(0)
        s = mp.mpc(1)
        previous_argument: mp.mpf | None = None
        last: tuple[mp.mpc, mp.mpc, mp.mpc, mp.mpc] | None = None
        for depth in range(1, 500):
            w_next = mp.log(w) / self.c
            v_next = v / (w * self.c) - w_next
            s_next = s / (w * self.c)
            w, v, s = w_next, v_next, s_next
            delta = w - self.L

            raw_argument = mp.arg(delta)
            if previous_argument is None:
                argument = raw_argument
            else:
                predicted = previous_argument - mp.arg(self.lam)
                turns = mp.nint((predicted - raw_argument) / (2 * mp.pi))
                argument = raw_argument + 2 * mp.pi * turns
            previous_argument = argument

            log_s = (
                depth * self.log_lam
                + mp.log(abs(delta))
                + 1j * argument
            )
            dlog_s = depth * self.dlog_lam + (v - self.dL) / delta
            slog_derivative = s / delta
            chi = log_s / self.log_lam
            gamma = (dlog_s - chi * self.dlog_lam) / slog_derivative
            last = (gamma / y_prime, log_s, dlog_s, slog_derivative)
            if abs(delta) < KOENIGS_STOP:
                assert last is not None
                field, log_s, dlog_s, slog_derivative = last
                return field, {
                    "depth": depth,
                    "terminal_distance": decimal(abs(delta)),
                    "log_S": complex_record(log_s),
                    "d_beta_log_S": complex_record(dlog_s),
                    "Sprime_over_S": complex_record(slog_derivative),
                }
        raise RuntimeError("Koenigs recursion did not reach the stopping radius")


def tail_sum(y: mp.mpf, derivative: mp.mpf, c: mp.mpf) -> tuple[mp.mpf, dict]:
    total = y / derivative
    current_y = y
    current_derivative = derivative
    last_term = abs(total)
    terms = 1
    for _ in range(1, 32):
        term = 1 / (c * current_derivative)
        total += term
        terms += 1
        last_term = abs(term)
        if last_term < mp.mpf("1e-65"):
            break
        current_y = mp.exp(c * current_y)
        current_derivative = c * current_y * current_derivative
    return total, {
        "terms": terms,
        "last_term": decimal(last_term),
        "remainder_bound": decimal(2 * last_term),
    }


def nested_exponential_expression(y: mp.mpf, c: mp.mpf, depth: int) -> str:
    expression = f"({decimal(y, 68)})"
    c_text = decimal(c, 68)
    for _ in range(depth):
        expression = f"exp(({c_text})*({expression}))"
    return expression


def raw_inverse_expressions(
    state_values: list[mp.mpf],
    anchor_c: mp.mpf,
    target_c: mp.mpf,
) -> tuple[list[str], mp.mpf, mp.mpf]:
    alpha = anchor_c / target_c
    offset = mp.log(alpha) / target_c
    expressions = []
    for value in state_values:
        tower = nested_exponential_expression(value, target_c, RAW_DEPTH)
        expressions.append(
            f"slog((({tower})-({decimal(offset, 68)}))"
            f"/({decimal(alpha, 68)}))-({RAW_DEPTH})"
        )
    return expressions, alpha, offset


def engine_state(
    engine: FatouGP,
    base_text: str,
    x_nodes: list[mp.mpf],
) -> tuple[list[mp.mpf], list[mp.mpf], list[mp.mpf], mp.mpf]:
    expressions: list[str] = []
    for x in x_nodes:
        text = decimal(x, 68)
        expressions.extend([f"sexp({text})", f"derivnum(t={text},sexp(t))"])
    values = engine.eval_batch(base_text, expressions)
    state = [mp.re(values[2 * j]) for j in range(len(x_nodes))]
    state_derivative = [mp.re(values[2 * j + 1]) for j in range(len(x_nodes))]
    coefficients = cheb_coefficients(state)
    derivative_coefficients = cheb_derivative(coefficients)
    direct_error = max(
        abs(cheb_eval(derivative_coefficients, x_nodes[j]) - state_derivative[j])
        for j in range(len(x_nodes))
    )
    return state, coefficients, derivative_coefficients, direct_error


def high_accuracy_rh(
    beta: mp.mpf,
    state_coefficients: list[mp.mpf],
    state_derivative_coefficients: list[mp.mpf],
    quadrature_nodes: list[mp.mpf],
    quadrature_weights: list[mp.mpf],
    target_nodes: list[mp.mpf],
) -> dict:
    koenigs = KoenigsData(beta)
    quadrature_u: list[mp.mpf] = []
    depth_min = 10**9
    depth_max = 0
    terminal_max = mp.mpf(0)
    cocycle_samples = []

    for index, x in enumerate(quadrature_nodes):
        y = mp.re(cheb_eval(state_coefficients, x))
        y_prime = mp.re(cheb_eval(state_derivative_coefficients, x))
        field, trace = koenigs.geometric_field(y, y_prime)
        quadrature_u.append(mp.im(field))
        depth_min = min(depth_min, trace["depth"])
        depth_max = max(depth_max, trace["depth"])
        terminal_max = max(terminal_max, mp.mpf(trace["terminal_distance"]))
        if index in (0, len(quadrature_nodes) // 2, len(quadrature_nodes) - 1):
            cocycle_samples.append(trace)

    geometric: list[mp.mpc] = []
    target_u: list[mp.mpf] = []
    for x in target_nodes:
        if x == 0:
            field = mp.mpc(0)
        elif x == 1:
            y = mp.re(cheb_eval(state_coefficients, x))
            y_prime = mp.re(cheb_eval(state_derivative_coefficients, x))
            field = mp.mpc(-koenigs.base * koenigs.c / y_prime)
        else:
            y = mp.re(cheb_eval(state_coefficients, x))
            y_prime = mp.re(cheb_eval(state_derivative_coefficients, x))
            field, _ = koenigs.geometric_field(y, y_prime)
        geometric.append(field)
        target_u.append(mp.im(field))

    hilbert: list[mp.mpf] = []
    for x, u_x in zip(target_nodes, target_u):
        value = mp.mpf(0)
        for t, weight, u_t in zip(
            quadrature_nodes, quadrature_weights, quadrature_u
        ):
            value += weight * (u_t - u_x) * mp.cot(mp.pi * (x - t))
        hilbert.append(value)

    unnormalized = [
        mp.re(geometric[j]) + hilbert[j] for j in range(len(target_nodes))
    ]
    normalization = unnormalized[target_nodes.index(mp.mpf(0))]
    selected = [value - normalization for value in unnormalized]

    tails = []
    kernel = []
    for x, selected_value in zip(target_nodes, selected):
        y = mp.re(cheb_eval(state_coefficients, x))
        y_prime = mp.re(cheb_eval(state_derivative_coefficients, x))
        tail, tail_trace = tail_sum(y, y_prime, koenigs.c)
        tails.append({"value": decimal(tail), **tail_trace})
        kernel.append(selected_value - tail)

    return {
        "selected_values": selected,
        "kernel_values": kernel,
        "geometric_values": geometric,
        "hilbert_values": hilbert,
        "tail_values": tails,
        "normalization_constant": normalization,
        "koenigs": {
            "L": complex_record(koenigs.L),
            "lambda": complex_record(koenigs.lam),
            "d_beta_L": complex_record(koenigs.dL),
            "d_beta_lambda": complex_record(koenigs.dlam),
            "depth_min": depth_min,
            "depth_max": depth_max,
            "terminal_distance_max": decimal(terminal_max),
            "sample_traces": cocycle_samples,
        },
    }


def stencil_diagnostic(
    beta: mp.mpf,
    x_nodes: list[mp.mpf],
    state_derivative_coefficients: list[mp.mpf],
    rh_selected: list[mp.mpf],
    gp_exe: Path,
    fatou_gp: str,
) -> dict:
    plus_beta = beta + STENCIL_STEP
    minus_beta = beta - STENCIL_STEP
    plus_base = decimal(mp.exp(mp.exp(plus_beta)), 72)
    minus_base = decimal(mp.exp(mp.exp(minus_beta)), 72)
    plus = FatouGP(
        gp_exe=gp_exe,
        dps=DPS,
        fatou_gp=fatou_gp,
        state_cache=False,
        persistent=True,
        init_timeout=1200,
        eval_timeout=300,
    )
    minus = FatouGP(
        gp_exe=gp_exe,
        dps=DPS,
        fatou_gp=fatou_gp,
        state_cache=False,
        persistent=True,
        init_timeout=1200,
        eval_timeout=300,
    )
    # exp-075c: die beiden Engines wurden nie geschlossen. stencil_diagnostic
    # laeuft einmal pro Knoten, bei sieben Knoten blieben also 14 gp.exe-
    # Prozesse bis zum Ende der Sitzung stehen -- jeder mit eigenem PARI-Stack.
    # Sie lassen sich NICHT aus der Knotenschleife heben, weil ihre Basen
    # exp(exp(beta +- STENCIL_STEP)) vom Knoten abhaengen; was fehlte, war das
    # Aufraeumen.
    expressions = [f"sexp({decimal(x, 68)})" for x in x_nodes]
    try:
        plus_values = plus.eval_batch(plus_base, expressions)
        minus_values = minus.eval_batch(minus_base, expressions)
    finally:
        for engine in (plus, minus):
            try:
                engine.close()
            except Exception:  # Aufraeumen darf die Diagnose nicht versenken
                pass
    direct = []
    differences = []
    for index, x in enumerate(x_nodes):
        derivative = mp.re(cheb_eval(state_derivative_coefficients, x))
        velocity = -mp.re(
            (plus_values[index] - minus_values[index]) / (2 * STENCIL_STEP)
        ) / derivative
        direct.append(velocity)
        differences.append(abs(velocity - rh_selected[index]))
    return {
        "step": decimal(STENCIL_STEP),
        "direct_values": [decimal(value) for value in direct],
        "max_abs_difference": decimal(max(differences)),
        "target_engine_used_only_for_diagnostic": True,
    }


def source_manifest(paths: Iterable[Path], root: Path) -> list[dict]:
    records = []
    for path in paths:
        records.append(
            {
                "path": path.relative_to(root).as_posix(),
                "sha256": sha256(path),
                "size_bytes": path.stat().st_size,
            }
        )
    return records


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gp-exe", type=Path, required=True)
    parser.add_argument("--fatou-gp", default="fork")
    args = parser.parse_args()

    mp.mp.dps = DPS
    repo = Path(__file__).resolve().parents[2]
    args.output.parent.mkdir(parents=True, exist_ok=True)

    beta_target = mp.log(mp.log(2))
    beta_a = mp.mpf(7) * beta_target / 8
    beta_b = mp.mpf(11) * beta_target / 12
    base_coordinate_nodes = [
        (1 - mp.cos(mp.pi * j / (BASE_N - 1))) / 2 for j in range(BASE_N)
    ]
    beta_nodes = [
        beta_a + (beta_b - beta_a) * tau for tau in base_coordinate_nodes
    ]
    base_nodes = [mp.exp(mp.exp(beta)) for beta in beta_nodes]
    x_state_nodes = cheb_lobatto_nodes(STATE_N)
    x_uniform_nodes = [mp.mpf(j) / RH_N for j in range(RH_N)]
    x_rh_targets = list(dict.fromkeys(x_uniform_nodes + x_state_nodes))
    quadrature_nodes, quadrature_weights, panels = panel_rule()

    engine = FatouGP(
        gp_exe=args.gp_exe,
        dps=DPS,
        fatou_gp=args.fatou_gp,
        state_cache=False,
        persistent=True,
        init_timeout=1200,
        eval_timeout=300,
    )
    anchor_engine = FatouGP(
        gp_exe=args.gp_exe,
        dps=DPS,
        fatou_gp=args.fatou_gp,
        state_cache=False,
        persistent=True,
        init_timeout=1200,
        eval_timeout=300,
    )
    anchor_base_text = decimal(base_nodes[0], 72)
    anchor_c = mp.exp(beta_nodes[0])

    node_records = []
    h_uniform_rows: list[list[mp.mpf]] = []
    psi_uniform_rows: list[list[mp.mpf]] = []

    for node_index, (tau, beta, base) in enumerate(
        zip(base_coordinate_nodes, beta_nodes, base_nodes)
    ):
        print(
            f"[8r1] node {node_index + 1}/{BASE_N}: "
            f"beta={decimal(beta, 18)} base={decimal(base, 18)}",
            flush=True,
        )
        base_text = decimal(base, 72)
        state_values, state_coeffs, state_derivative_coeffs, derivative_gate = (
            engine_state(engine, base_text, x_state_nodes)
        )

        inverse_expressions, alpha, offset = raw_inverse_expressions(
            state_values, anchor_c, mp.exp(beta)
        )
        inverse_values = anchor_engine.eval_batch(
            anchor_base_text, inverse_expressions
        )
        U_values = [mp.re(value) for value in inverse_values]
        zero_index = x_state_nodes.index(mp.mpf(0))
        q_value = U_values[zero_index]
        h_values = [
            U_values[j] - x_state_nodes[j] - q_value
            for j in range(STATE_N)
        ]
        h_coeffs = cheb_coefficients(h_values)
        h_derivative_coeffs = cheb_derivative(h_coeffs)
        h_uniform = [
            mp.re(cheb_eval(h_coeffs, x)) for x in x_uniform_nodes
        ]
        h_prime_uniform = [
            mp.re(cheb_eval(h_derivative_coeffs, x)) for x in x_uniform_nodes
        ]

        rh = high_accuracy_rh(
            beta,
            state_coeffs,
            state_derivative_coeffs,
            quadrature_nodes,
            quadrature_weights,
            x_rh_targets,
        )
        target_index = {value: index for index, value in enumerate(x_rh_targets)}
        selected_uniform = [
            rh["selected_values"][target_index[x]] for x in x_uniform_nodes
        ]
        kernel_uniform = [
            rh["kernel_values"][target_index[x]] for x in x_uniform_nodes
        ]
        selected_state = [
            rh["selected_values"][target_index[x]] for x in x_state_nodes
        ]
        kernel_state = [
            rh["kernel_values"][target_index[x]] for x in x_state_nodes
        ]
        kernel_modes = fourier_coefficients(kernel_uniform)
        selected_modes = fourier_coefficients(selected_uniform)

        q_prime = -(1 + h_prime_uniform[0]) * kernel_uniform[0]
        psi_uniform = [
            -(1 + h_prime_uniform[j]) * kernel_uniform[j] - q_prime
            for j in range(RH_N)
        ]

        diagnostic = stencil_diagnostic(
            beta,
            x_uniform_nodes,
            state_derivative_coeffs,
            selected_uniform,
            args.gp_exe,
            args.fatou_gp,
        )
        print(
            f"[8r1] node {node_index + 1}: "
            f"RH/stencil={diagnostic['max_abs_difference']}",
            flush=True,
        )

        state_fourier = fourier_coefficients(h_uniform)
        tail_window = [
            abs(kernel_modes[k]) for k in range(RH_HEAD + 1, RH_N // 2)
        ]
        measured_tail = 2 * sum(tail_window)
        last_ratio = max(
            (
                abs(kernel_modes[k]) / abs(kernel_modes[k - 1])
                for k in range(RH_HEAD + 2, RH_N // 2)
                if abs(kernel_modes[k - 1]) != 0
            ),
            default=mp.mpf("0.5"),
        )
        continuation_ratio = min(mp.mpf("0.95"), 2 * last_ratio)
        continuation = (
            2
            * abs(kernel_modes[RH_N // 2 - 1])
            * continuation_ratio
            / (1 - continuation_ratio)
        )
        fourier_tail_bound = measured_tail + continuation
        cheb_tail_bound = 4 * sum(abs(value) for value in state_coeffs[-6:])

        node_records.append(
            {
                "index": node_index,
                "tau": decimal(tau),
                "beta": decimal(beta),
                "base": decimal(base),
                "raw_transform": {
                    "anchor_base": decimal(base_nodes[0]),
                    "anchor_log_base": decimal(anchor_c),
                    "target_log_base": decimal(mp.exp(beta)),
                    "alpha": decimal(alpha),
                    "offset": decimal(offset),
                    "ladder_depth": RAW_DEPTH,
                    "q": decimal(q_value),
                    "U_values_at_state_nodes": [
                        decimal(value) for value in U_values
                    ],
                },
                "state": {
                    "representation": "Chebyshev on x in [0,1], s=2*x-1",
                    "values": [decimal(value) for value in state_values],
                    "coefficients": [decimal(value) for value in state_coeffs],
                    "derivative_coefficients": [
                        decimal(value) for value in state_derivative_coeffs
                    ],
                    "derivative_interpolation_gate": decimal(derivative_gate),
                    "chebyshev_tail_bound": decimal(cheb_tail_bound),
                },
                "normalized_hub": {
                    "h_values_at_state_nodes": [
                        decimal(value) for value in h_values
                    ],
                    "h_coefficients": [decimal(value) for value in h_coeffs],
                    "h_derivative_coefficients": [
                        decimal(value) for value in h_derivative_coeffs
                    ],
                    "h_uniform_values": [
                        decimal(value) for value in h_uniform
                    ],
                    "h_fourier_modes": [
                        {
                            "k": k,
                            **complex_record(state_fourier[k]),
                        }
                        for k in range(-RH_HEAD, RH_HEAD + 1)
                    ],
                    "q_prime_from_normalization": decimal(q_prime),
                },
                "rh": {
                    "selected_values_uniform": [
                        decimal(value) for value in selected_uniform
                    ],
                    "kernel_values_uniform": [
                        decimal(value) for value in kernel_uniform
                    ],
                    "selected_values_state_nodes": [
                        decimal(value) for value in selected_state
                    ],
                    "kernel_values_state_nodes": [
                        decimal(value) for value in kernel_state
                    ],
                    "selected_modes": [
                        {"k": k, **complex_record(selected_modes[k])}
                        for k in range(-RH_HEAD, RH_HEAD + 1)
                    ],
                    "kernel_modes": [
                        {"k": k, **complex_record(kernel_modes[k])}
                        for k in range(-RH_HEAD, RH_HEAD + 1)
                    ],
                    "analytic_tail": {
                        "measured_resolved_tail": decimal(measured_tail),
                        "continuation_ratio": decimal(continuation_ratio),
                        "continuation_bound": decimal(continuation),
                        "total_bound": decimal(fourier_tail_bound),
                    },
                    "normalization_constant": decimal(
                        rh["normalization_constant"]
                    ),
                    "koenigs": rh["koenigs"],
                    "tail_values": [
                        rh["tail_values"][target_index[x]]
                        for x in x_uniform_nodes
                    ],
                    "engine_stencil_diagnostic": diagnostic,
                },
            }
        )
        h_uniform_rows.append(h_uniform)
        psi_uniform_rows.append(psi_uniform)

    # Differentiate the seven-node base path independently of the RH values.
    base_s_nodes = [2 * tau - 1 for tau in base_coordinate_nodes]
    base_diff = differentiation_matrix(base_s_nodes)
    path_derivative = [
        [mp.mpf(0) for _ in range(RH_N)] for _ in range(BASE_N)
    ]
    for i in range(BASE_N):
        for x_index in range(RH_N):
            derivative_s = sum(
                base_diff[i][j] * h_uniform_rows[j][x_index]
                for j in range(BASE_N)
            )
            path_derivative[i][x_index] = (
                2 * derivative_s / (beta_b - beta_a)
            )

    residual_rows = [
        [
            path_derivative[i][j] - psi_uniform_rows[i][j]
            for j in range(RH_N)
        ]
        for i in range(BASE_N)
    ]
    residual_modes = [fourier_coefficients(row) for row in residual_rows]

    gp_path = Path(args.gp_exe).resolve()
    fatou_path = (
        repo / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
        if args.fatou_gp == "fork"
        else Path(args.fatou_gp).resolve()
    )
    producer_path = Path(__file__).resolve()
    convention_blob = {
        "beta": "log(log(base))",
        "segment": {
            "beta_a": "7/8*log(log(2))",
            "beta_b": "11/12*log(log(2))",
        },
        "state_chebyshev_count": STATE_N,
        "base_lobatto_count": BASE_N,
        "rh_uniform_count": RH_N,
        "rh_head": RH_HEAD,
        "panel_order": PANEL_ORDER,
        "panel_grading": "decades, symmetric at 1/2",
        "tmin": decimal(TMIN),
        "koenigs_stop": decimal(KOENIGS_STOP),
        "raw_depth": RAW_DEPTH,
    }
    convention_sha = hashlib.sha256(
        json.dumps(convention_blob, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    payload = {
        "schema_version": "paper-vi-trajectory-8r1-v1",
        "segment_id": "8r1",
        "status": "recomputed-candidate",
        "provenance": {
            "original_files_found": False,
            "candidate_target_engine_used": True,
            "rh_operator_target_engine_used": False,
            "target_engine_use": (
                "state centers and blind beta-stencil diagnostic only; "
                "the RH field is rebuilt from Lambert-W/Koenigs data"
            ),
            "lineage": (
                "direct Kneser centers -> raw re-anchor transform -> "
                "normalized hub; independent Koenigs RH reconstruction"
            ),
        },
        "segment": {
            "beta_a": decimal(beta_a),
            "beta_b": decimal(beta_b),
            "base_a": decimal(base_nodes[0]),
            "base_b": decimal(base_nodes[-1]),
            "absolute_beta_length": decimal(abs(beta_b - beta_a)),
        },
        "precision": {
            "decimal_digits": DPS,
            "pari_gp": subprocess.check_output(
                [str(gp_path), "-f", "-q"],
                input=b"print(version())\n\\q\n",
            )
            .decode(errors="replace")
            .strip(),
            "mpmath": mp.__version__,
            "python": platform.python_version(),
        },
        "conventions": convention_blob,
        "conventions_sha256": convention_sha,
        "transforms": {
            "state_chebyshev_nodes": [
                decimal(value) for value in x_state_nodes
            ],
            "rh_uniform_nodes": [decimal(value) for value in x_uniform_nodes],
            "base_lobatto_tau_nodes": [
                decimal(value) for value in base_coordinate_nodes
            ],
            "base_lobatto_beta_nodes": [decimal(value) for value in beta_nodes],
            "base_barycentric_weights": [
                decimal(value) for value in barycentric_weights(base_s_nodes)
            ],
            "base_differentiation_matrix_ds": [
                [decimal(value) for value in row] for row in base_diff
            ],
            "quadrature": {
                "panels": [
                    {"a": decimal(a), "b": decimal(b)} for a, b in panels
                ],
                "nodes": [decimal(value) for value in quadrature_nodes],
                "weights": [decimal(value) for value in quadrature_weights],
            },
        },
        "nodes": node_records,
        "volterra": {
            "path_derivative_uniform": [
                [decimal(value) for value in row] for row in path_derivative
            ],
            "rhs_uniform": [
                [decimal(value) for value in row] for row in psi_uniform_rows
            ],
            "residual_uniform": [
                [decimal(value) for value in row] for row in residual_rows
            ],
            "residual_modes": [
                [
                    {"k": k, **complex_record(modes[k])}
                    for k in range(-RH_HEAD, RH_HEAD + 1)
                ]
                for modes in residual_modes
            ],
            "max_point_residual": decimal(
                max(abs(value) for row in residual_rows for value in row)
            ),
        },
        "source_manifest": source_manifest(
            [producer_path, fatou_path, gp_path], repo
            if gp_path.is_relative_to(repo)
            else Path("/"),
        )
        if gp_path.is_relative_to(repo)
        # exp-075: .as_posix() auch in diesem Zweig. Er greift, sobald die
        # gp-Binaerdatei ausserhalb des Repos liegt -- unter Windows also im
        # Normalfall -- und schrieb dann plattformeigene Trennzeichen ins
        # Manifest. Zweck des Manifests ist der Vergleich ueber Maschinen
        # hinweg; zwei Eintraege, die sich nur in / gegen \ unterscheiden,
        # sehen dabei aus wie verschiedene Dateien.
        else [
            {
                "path": producer_path.relative_to(repo).as_posix(),
                "sha256": sha256(producer_path),
                "size_bytes": producer_path.stat().st_size,
            },
            {
                "path": fatou_path.relative_to(repo).as_posix(),
                "sha256": sha256(fatou_path),
                "size_bytes": fatou_path.stat().st_size,
            },
            {
                "path": gp_path.as_posix(),
                "sha256": sha256(gp_path),
                "size_bytes": gp_path.stat().st_size,
            },
        ],
    }

    # Hash the canonical content without the self-hash, then add it.
    canonical = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    payload["payload_sha256"] = hashlib.sha256(canonical).hexdigest()
    args.output.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"[8r1] wrote {args.output}", flush=True)
    print(f"[8r1] payload sha256 {payload['payload_sha256']}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
