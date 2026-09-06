"""Outward-rounded interval producer for the Paper-VI segment ``8r1``.

This program does not consume the floating-point RH values stored in the
trajectory as proof input.  It rebuilds the primitive Koenigs, quotient and
tail arrays with Arb balls.  Dependence on the base parameter is retained as
an ``acb_series`` Taylor model on each consecutive Lobatto panel, so that the
independent replay can perform the DFT, Hardy projection, normalization,
Volterra residual and Lipschitz majorants *before* evaluating the base box.

The producer and replay intentionally share no Python module.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable

from flint import acb, acb_poly, acb_series, arb, ctx


ARBITRARY_PRECISION_BITS = 448
BETA_SERIES_ORDER = 22
STATE_COUNT = 48
RH_COUNT = 96
HUB_MODE_CUTOFF = 16
KOENIGS_STOP = arb("1e-16")
KOENIGS_MAX_DEPTH = 360
KOENIGS_GAMMA_FACTOR = arb(128)
KOENIGS_GAMMA_PRIME_FACTOR = arb(1024)
BETA_SERIES_REMAINDER = arb("2e-16")
TAIL_STOP = arb("1e-50")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def finite_ball(value: arb) -> bool:
    return bool(value.is_finite())


def outward_float_bounds(value: arb) -> list[str]:
    """Return binary64-facing decimal endpoints, widened away from the ball."""

    if not finite_ball(value):
        raise RuntimeError(f"non-finite Arb ball: {value}")
    lower = math.nextafter(float(value.lower()), -math.inf)
    upper = math.nextafter(float(value.upper()), math.inf)
    return [format(lower, ".17e"), format(upper, ".17e")]


def encode_acb(value: acb) -> dict[str, list[str]]:
    return {
        "re": outward_float_bounds(value.real),
        "im": outward_float_bounds(value.imag),
    }


def encode_series(value: acb_series) -> list[dict[str, list[str]]]:
    coefficients = value.coeffs()
    zero = acb(0)
    return [
        encode_acb(coefficients[index] if index < len(coefficients) else zero)
        for index in range(BETA_SERIES_ORDER)
    ]


def series_box(value: acb_series, radius: str = "1") -> acb:
    return acb_poly(value.coeffs())(arb(f"0 +/- {radius}"))


def add_complex_error(value: acb_series, error: arb) -> acb_series:
    coefficients = value.coeffs()
    while len(coefficients) < BETA_SERIES_ORDER:
        coefficients.append(acb(0))
    coefficients[0] += acb(arb(0, error), arb(0, error))
    return acb_series(coefficients, BETA_SERIES_ORDER)


def cheb_coefficients(values: list[arb | acb]) -> list[arb | acb]:
    n = len(values) - 1
    coefficients: list[arb | acb] = []
    for k in range(n + 1):
        total = values[0] / 2 + values[-1] * ((-1) ** k) / 2
        for j in range(1, n):
            total += values[j] * (arb.pi() * j * k / n).cos()
        coefficient = 2 * total / n
        if k in (0, n):
            coefficient /= 2
        coefficients.append(coefficient)
    return coefficients


def cheb_eval(coefficients, point):
    zero = (
        acb_series([], BETA_SERIES_ORDER)
        if isinstance(point, acb_series)
        or any(isinstance(value, acb_series) for value in coefficients)
        else acb(0)
    )
    b1 = zero
    b2 = zero
    for coefficient in coefficients[:0:-1]:
        b0 = 2 * point * b1 - b2 + coefficient
        b2, b1 = b1, b0
    return point * b1 - b2 + coefficients[0]


def state_second_derivative(
    derivative_coefficients: list[arb],
) -> list[arb]:
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


def base_model(
    values_at_nodes: Iterable[arb | acb],
    base_series: acb_series,
) -> acb_series:
    # Stored base nodes are ascending -cos(j*pi/6); DCT order is descending.
    coefficients = cheb_coefficients(list(reversed(list(values_at_nodes))))
    return cheb_eval(coefficients, base_series)


def differentiate_beta_series(
    value: acb_series, beta_halfwidth: arb
) -> acb_series:
    coefficients = value.coeffs()
    derivative = [
        (index + 1) * coefficients[index + 1] / beta_halfwidth
        for index in range(len(coefficients) - 1)
    ]
    return acb_series(derivative, BETA_SERIES_ORDER)


def dft(values: list[acb]) -> dict[int, acb]:
    count = len(values)
    output: dict[int, acb] = {}
    for k in range(-count // 2 + 1, count // 2):
        total = acb(0)
        for j, value in enumerate(values):
            angle = -2 * arb.pi() * k * j / count
            total += value * acb(angle.cos(), angle.sin())
        output[k] = total / count
    return output


def build_hub_modes(trajectory: dict) -> dict[int, list[acb]]:
    """Return exact outward DFT heads at the seven base nodes."""

    result = {
        k: [] for k in range(-HUB_MODE_CUTOFF, HUB_MODE_CUTOFF + 1) if k
    }
    for node in trajectory["nodes"]:
        values = [acb(arb(text)) for text in node["normalized_hub"]["h_uniform_values"]]
        modes = dft(values)
        for k in result:
            result[k].append(modes[k])
    return result


def koenigs_models(
    beta: acb_series,
    y: acb_series,
    want_prime: bool = True,
) -> tuple[acb_series, acb_series | None, dict]:
    """Compute Gamma and, on request, dGamma/dy with a stopped recursion.

    exp-075d: der Hilbert-Zertifizierer verwirft dGamma/dy (`gamma, _, trace`
    in certify_hilbert_8r1.py), liess es aber vollstaendig ausrechnen -- samt
    der u/s2-Rekursion, die AUSSCHLIESSLICH dafuer existiert und pro Iteration
    zwei weitere Reihenoperationen kostet. Mit want_prime=False entfaellt sie.

    Der zurueckgegebene Trace bleibt dabei BITGLEICH: gamma_prime_error haengt
    nur an delta_box, denominator und den Konstanten, nicht an gamma_prime.
    """

    c = beta.exp()
    L = -(-c).lambertw(-1) / c
    lam = c * L
    log_lam = lam.log()
    dL = c * L * L / (1 - lam)
    dlam = c * (dL + L)
    dlog_lam = dlam / lam

    w = y
    v = acb_series([], BETA_SERIES_ORDER)
    s = acb_series([1], BETA_SERIES_ORDER)
    u = acb_series([], BETA_SERIES_ORDER) if want_prime else None
    s2 = acb_series([], BETA_SERIES_ORDER) if want_prime else None
    previous_argument: arb | None = None
    minimum_phase_margin: arb | None = None

    for depth in range(1, KOENIGS_MAX_DEPTH + 1):
        w_next = w.log() / c
        s_next = s / (c * w)
        v_next = v / (c * w) - w_next
        if want_prime:
            u_next = u / (c * w) - v * s / (c * w * w) - s_next
            s2_next = s2 / (c * w) - s * s / (c * w * w)
            w, v, s, u, s2 = w_next, v_next, s_next, u_next, s2_next
        else:
            w, v, s = w_next, v_next, s_next
        delta = w - L

        raw_argument = delta.coeffs()[0].arg()
        if previous_argument is None:
            turn = 0
            adjusted_argument = raw_argument
        else:
            predicted = previous_argument - lam.coeffs()[0].arg()
            turn = round(
                float(((predicted - raw_argument) / (2 * arb.pi())).mid())
            )
            adjusted_argument = raw_argument + 2 * arb.pi() * turn
            margin = arb.pi() - (predicted - adjusted_argument).abs_upper()
            if margin <= 0:
                raise RuntimeError("Koenigs phase turn is not unique")
            if minimum_phase_margin is None or margin < minimum_phase_margin:
                minimum_phase_margin = margin
        previous_argument = adjusted_argument

        # exp-075d: der Gamma/Gamma'-Block stand frueher hier, wurde also in
        # JEDER Iteration gerechnet -- bei den gemessenen ~178 Tiefen der
        # graded Knoten also 177-mal umsonst, denn benutzt wird er nur im
        # Abbruchzweig, der sofort zurueckkehrt.
        # Die Verschiebung ist bitgleich per Konstruktion: keine der zehn
        # Groessen (log_s, dlog_s, slog_derivative, chi, numerator, gamma,
        # d_dlog_s, d_slog_derivative, d_numerator, gamma_prime) wird in der
        # naechsten Iteration gelesen; schleifengetragen sind nur die
        # Rekursion selbst, previous_argument und minimum_phase_margin, und
        # die bleiben oben. In der Abbruchiteration liegen depth, delta, v,
        # dL, s, u, s2 und turn unveraendert vor.
        delta_box = series_box(delta)
        if delta_box.abs_upper() < KOENIGS_STOP:
            log_s = (
                depth * log_lam
                + delta.log()
                + acb(0, 2 * arb.pi() * turn)
            )
            dlog_s = depth * dlog_lam + (v - dL) / delta
            slog_derivative = s / delta
            chi = log_s / log_lam
            numerator = dlog_s - chi * dlog_lam
            gamma = numerator / slog_derivative

            if want_prime:
                d_dlog_s = (u * delta - (v - dL) * s) / (delta * delta)
                d_slog_derivative = (s2 * delta - s * s) / (delta * delta)
                d_numerator = (
                    d_dlog_s - (slog_derivative / log_lam) * dlog_lam
                )
                gamma_prime = (
                    d_numerator * slog_derivative
                    - numerator * d_slog_derivative
                ) / (slog_derivative * slog_derivative)
            else:
                gamma_prime = None

            inv_lambda = series_box(1 / lam).abs_upper()
            if inv_lambda >= 1:
                raise RuntimeError("local Koenigs contraction is not strict")
            denominator = 1 - inv_lambda
            gamma_error = (
                KOENIGS_GAMMA_FACTOR
                * delta_box.abs_upper()
                / denominator**3
                + BETA_SERIES_REMAINDER
            )
            gamma_prime_error = (
                KOENIGS_GAMMA_PRIME_FACTOR
                * delta_box.abs_upper()
                / denominator**5
                + BETA_SERIES_REMAINDER
            )
            gamma = add_complex_error(gamma, gamma_error)
            if gamma_prime is not None:
                gamma_prime = add_complex_error(gamma_prime, gamma_prime_error)
            return gamma, gamma_prime, {
                "depth": depth,
                "terminal_distance": outward_float_bounds(
                    delta_box.abs_upper()
                ),
                "inverse_lambda_upper": outward_float_bounds(inv_lambda),
                "phase_margin_lower": outward_float_bounds(
                    minimum_phase_margin or arb.pi()
                ),
                "gamma_remainder": outward_float_bounds(gamma_error),
                "gamma_prime_remainder": outward_float_bounds(
                    gamma_prime_error
                ),
            }
    raise RuntimeError("stopped Koenigs recursion did not contract")


def tail_models(
    c: acb_series,
    y: acb_series,
    y_prime: acb_series,
) -> tuple[acb_series, acb_series, acb_series, dict]:
    """Return S, partial_y S and partial_yprime S."""

    total = y / y_prime
    derivative_y = 1 / y_prime
    derivative_y_prime = -y / (y_prime * y_prime)

    current_y = y
    current_y_prime = y_prime
    dy_dy = acb_series([1], BETA_SERIES_ORDER)
    dyp_dy = acb_series([], BETA_SERIES_ORDER)
    dy_dyp = acb_series([], BETA_SERIES_ORDER)
    dyp_dyp = acb_series([1], BETA_SERIES_ORDER)
    last_term = total
    terms = 1

    for _ in range(1, 32):
        term = 1 / (c * current_y_prime)
        total += term
        derivative_y += -dyp_dy / (
            c * current_y_prime * current_y_prime
        )
        derivative_y_prime += -dyp_dyp / (
            c * current_y_prime * current_y_prime
        )
        last_term = term
        terms += 1
        if series_box(term).abs_upper() < TAIL_STOP:
            break

        next_y = (c * current_y).exp()
        next_dy_dy = c * next_y * dy_dy
        next_dy_dyp = c * next_y * dy_dyp
        next_y_prime = c * next_y * current_y_prime
        next_dyp_dy = c * (
            next_dy_dy * current_y_prime + next_y * dyp_dy
        )
        next_dyp_dyp = c * (
            next_dy_dyp * current_y_prime + next_y * dyp_dyp
        )
        (
            current_y,
            current_y_prime,
            dy_dy,
            dyp_dy,
            dy_dyp,
            dyp_dyp,
        ) = (
            next_y,
            next_y_prime,
            next_dy_dy,
            next_dyp_dy,
            next_dy_dyp,
            next_dyp_dyp,
        )
    else:
        # exp-075c: ohne diesen Zweig lief die Schleife bei ausbleibender
        # Kontraktion einfach aus und der Code rechnete unten mit
        # remainder = 2*|letzter Term| weiter. Der Faktor 2 unterstellt eine
        # geometrische Reihe mit Quotient unter 1/2 -- genau die Annahme, die
        # ohne den break NICHT nachgewiesen ist. Das Ergebnis waere ein
        # unbegruendetes Restglied ohne Fehlermeldung gewesen.
        # Die drei vergleichbaren Schleifen im Repo brechen an dieser Stelle
        # ab: certify_8r1.py:282, tails.py:172-173, koenigs.py:301.
        # Bei 8r1 greift der break nach 6-7 Termen, der Zweig ist dort also
        # unerreichbar und die Zahlen bleiben bitgleich.
        raise RuntimeError(
            "tail recursion did not reach TAIL_STOP within 32 terms; "
            "the factor-2 geometric remainder below would be unjustified"
        )

    remainder = 2 * series_box(last_term).abs_upper()
    total = add_complex_error(total, remainder + BETA_SERIES_REMAINDER)
    derivative_y = add_complex_error(
        derivative_y, 4 * remainder + BETA_SERIES_REMAINDER
    )
    derivative_y_prime = add_complex_error(
        derivative_y_prime, 4 * remainder + BETA_SERIES_REMAINDER
    )
    return total, derivative_y, derivative_y_prime, {
        "terms": terms,
        "last_term_upper": outward_float_bounds(
            series_box(last_term).abs_upper()
        ),
        "remainder": outward_float_bounds(remainder),
    }


def interval_hull(records: list[list[str]]) -> list[str]:
    lower = min(float(record[0]) for record in records)
    upper = max(float(record[1]) for record in records)
    return [format(lower, ".17e"), format(upper, ".17e")]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    ctx.prec = ARBITRARY_PRECISION_BITS

    # exp-075 (2026-07-31): python-flint kappt JEDE Serienoperation bei der globalen
    # Laenge ctx.cap, deren Vorgabe 10 ist. Keines der drei Programme hat sie je
    # gesetzt, waehrend die Payload "beta_series_order": 22 behauptete -- alle
    # beta-Modelle waren in Wahrheit Polynome vom Grad 9, die Koeffizienten 10..21
    # exakt null. Gemeldet im Closure-Dokument, auf dieser Maschine unabhaengig
    # nachgeprueft (flint 0.9.0, ctx.cap == 10). Ohne diese Zeile erzeugt der Lauf
    # ein Artefakt, dessen Etikett nicht stimmt.
    ctx.cap = BETA_SERIES_ORDER
    trajectory = json.loads(args.trajectory.read_text(encoding="utf-8"))
    if trajectory.get("schema_version") != "paper-vi-trajectory-8r1-v1":
        raise RuntimeError("unexpected trajectory schema")
    if len(trajectory.get("nodes", [])) != 7:
        raise RuntimeError("8r1 requires seven Lobatto states")

    beta_a = arb(trajectory["segment"]["beta_a"])
    beta_b = arb(trajectory["segment"]["beta_b"])
    beta_nodes = [
        arb(value)
        for value in trajectory["transforms"]["base_lobatto_beta_nodes"]
    ]
    hub_modes = build_hub_modes(trajectory)

    state_coefficients_at_nodes = [
        [arb(value) for value in node["state"]["coefficients"]]
        for node in trajectory["nodes"]
    ]
    state_derivatives_at_nodes = [
        [arb(value) for value in node["state"]["derivative_coefficients"]]
        for node in trajectory["nodes"]
    ]
    state_second_at_nodes = [
        state_second_derivative(row) for row in state_derivatives_at_nodes
    ]

    panels = []
    all_phase_margins: list[list[str]] = []
    all_inverse_lambda: list[list[str]] = []
    all_terminal_distances: list[list[str]] = []
    all_gamma_remainders: list[list[str]] = []
    all_gamma_prime_remainders: list[list[str]] = []

    for panel_index in range(6):
        beta_left = beta_nodes[panel_index]
        beta_right = beta_nodes[panel_index + 1]
        beta_midpoint = (beta_left + beta_right) / 2
        beta_halfwidth = (beta_right - beta_left) / 2
        beta_series = acb_series(
            [acb(beta_midpoint), acb(beta_halfwidth)],
            BETA_SERIES_ORDER,
        )
        base_series = (
            2 * (beta_series - beta_a) / (beta_b - beta_a) - 1
        )

        print(
            f"[8r1-interval] panel {panel_index + 1}/6 "
            f"beta={beta_left}..{beta_right}",
            flush=True,
        )

        state_models = []
        state_derivative_models = []
        state_second_models = []
        for coefficient_index in range(STATE_COUNT):
            state_models.append(
                base_model(
                    [
                        state_coefficients_at_nodes[node][coefficient_index]
                        for node in range(7)
                    ],
                    base_series,
                )
            )
            state_derivative_models.append(
                base_model(
                    [
                        state_derivatives_at_nodes[node][coefficient_index]
                        for node in range(7)
                    ],
                    base_series,
                )
            )
            state_second_models.append(
                base_model(
                    [
                        state_second_at_nodes[node][coefficient_index]
                        for node in range(7)
                    ],
                    base_series,
                )
            )

        hub_mode_models = {
            k: base_model(values, base_series)
            for k, values in hub_modes.items()
        }

        primitive_rows = []
        panel_phase_margins: list[list[str]] = []
        panel_inverse_lambda: list[list[str]] = []
        panel_terminal: list[list[str]] = []
        panel_gamma_remainder: list[list[str]] = []
        panel_gamma_prime_remainder: list[list[str]] = []

        for x_index in range(RH_COUNT):
            x = arb(x_index) / RH_COUNT
            state_s = acb(2 * x - 1)
            y = cheb_eval(state_models, state_s)
            y_prime = cheb_eval(state_derivative_models, state_s)
            y_second = cheb_eval(state_second_models, state_s)

            h = acb_series([], BETA_SERIES_ORDER)
            h_prime = acb_series([], BETA_SERIES_ORDER)
            h_second = acb_series([], BETA_SERIES_ORDER)
            for k, mode in hub_mode_models.items():
                angle = 2 * arb.pi() * k * x
                exponential = acb(angle.cos(), angle.sin())
                h += mode * (exponential - 1)
                frequency = acb(0, 2 * arb.pi() * k)
                h_prime += mode * frequency * exponential
                h_second += (
                    mode * (-(2 * arb.pi() * k) ** 2) * exponential
                )
            h_beta = differentiate_beta_series(h, beta_halfwidth)

            if x_index == 0:
                gamma = acb_series([], BETA_SERIES_ORDER)
                gamma_prime = acb_series([], BETA_SERIES_ORDER)
                koenigs_trace = {
                    "depth": 0,
                    "terminal_distance": outward_float_bounds(arb(0)),
                    "inverse_lambda_upper": outward_float_bounds(arb(0)),
                    "phase_margin_lower": outward_float_bounds(arb.pi()),
                    "gamma_remainder": outward_float_bounds(arb(0)),
                    "gamma_prime_remainder": outward_float_bounds(arb(0)),
                    "endpoint_stopped_model": True,
                }
            else:
                gamma, gamma_prime, koenigs_trace = koenigs_models(
                    beta_series, y
                )

            tail, tail_y, tail_y_prime, tail_trace = tail_models(
                beta_series.exp(), y, y_prime
            )
            primitive_rows.append(
                {
                    "x_index": x_index,
                    "T": encode_series(y),
                    "T_prime": encode_series(y_prime),
                    "T_second": encode_series(y_second),
                    "h": encode_series(h),
                    "h_prime": encode_series(h_prime),
                    "h_second": encode_series(h_second),
                    "h_beta": encode_series(h_beta),
                    "Gamma": encode_series(gamma),
                    "Gamma_prime": encode_series(gamma_prime),
                    "tail": encode_series(tail),
                    "tail_partial_T": encode_series(tail_y),
                    "tail_partial_T_prime": encode_series(tail_y_prime),
                    "koenigs_trace": koenigs_trace,
                    "tail_trace": tail_trace,
                }
            )
            panel_phase_margins.append(
                koenigs_trace["phase_margin_lower"]
            )
            panel_inverse_lambda.append(
                koenigs_trace["inverse_lambda_upper"]
            )
            panel_terminal.append(koenigs_trace["terminal_distance"])
            panel_gamma_remainder.append(
                koenigs_trace["gamma_remainder"]
            )
            panel_gamma_prime_remainder.append(
                koenigs_trace["gamma_prime_remainder"]
            )
            if x_index and x_index % 24 == 0:
                print(
                    f"[8r1-interval] panel {panel_index + 1}: "
                    f"x {x_index}/{RH_COUNT}",
                    flush=True,
                )

        panels.append(
            {
                "index": panel_index,
                "beta_left": outward_float_bounds(beta_left),
                "beta_right": outward_float_bounds(beta_right),
                "beta_midpoint": outward_float_bounds(beta_midpoint),
                "beta_halfwidth": outward_float_bounds(beta_halfwidth),
                "primitive_rows": primitive_rows,
            }
        )
        all_phase_margins.extend(panel_phase_margins)
        all_inverse_lambda.extend(panel_inverse_lambda)
        all_terminal_distances.extend(panel_terminal)
        all_gamma_remainders.extend(panel_gamma_remainder)
        all_gamma_prime_remainders.extend(panel_gamma_prime_remainder)

    producer_path = Path(__file__).resolve()
    payload = {
        "schema_version": "paper-vi-8r1-interval-primitives-v1",
        "segment_id": "8r1",
        "interval_backend": {
            "name": "python-flint/Arb",
            "version": __import__("flint").__version__,
            "precision_bits": ARBITRARY_PRECISION_BITS,
            "directed_rounding": True,
            "serialized_endpoints": "outward-widened IEEE-754 binary64",
        },
        "trajectory": {
            "path": args.trajectory.name,
            "file_sha256": sha256(args.trajectory),
            "payload_sha256": trajectory["payload_sha256"],
        },
        "conventions_sha256": trajectory["conventions_sha256"],
        "dimensions": {
            "base_lobatto_nodes": 7,
            "base_panels": 6,
            "beta_series_order": BETA_SERIES_ORDER,
            "state_chebyshev_count": STATE_COUNT,
            "rh_uniform_count": RH_COUNT,
            "hub_mode_cutoff": HUB_MODE_CUTOFF,
        },
        "proof_constants": {
            "koenigs_stop": outward_float_bounds(KOENIGS_STOP),
            "koenigs_gamma_factor": outward_float_bounds(
                KOENIGS_GAMMA_FACTOR
            ),
            "koenigs_gamma_prime_factor": outward_float_bounds(
                KOENIGS_GAMMA_PRIME_FACTOR
            ),
            "beta_series_remainder": outward_float_bounds(
                BETA_SERIES_REMAINDER
            ),
            "tail_stop": outward_float_bounds(TAIL_STOP),
            "endpoint_model": (
                "Gamma(T(0))=0; the removable stopped-coordinate endpoint "
                "model is used before any quotient"
            ),
        },
        "global_koenigs_gates": {
            "phase_margin_hull": interval_hull(all_phase_margins),
            "inverse_lambda_hull": interval_hull(all_inverse_lambda),
            "terminal_distance_hull": interval_hull(
                all_terminal_distances
            ),
            "gamma_remainder_hull": interval_hull(
                all_gamma_remainders
            ),
            "gamma_prime_remainder_hull": interval_hull(
                all_gamma_prime_remainders
            ),
        },
        "panels": panels,
        "source_manifest": [
            {
                "path": producer_path.name,
                "sha256": sha256(producer_path),
                "size_bytes": producer_path.stat().st_size,
            },
            {
                "path": args.trajectory.name,
                "sha256": sha256(args.trajectory),
                "size_bytes": args.trajectory.stat().st_size,
            },
        ],
    }
    canonical = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    payload["payload_sha256"] = hashlib.sha256(canonical).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, separators=(",", ":"), ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"[8r1-interval] wrote {args.output}", flush=True)
    print(
        f"[8r1-interval] payload sha256 {payload['payload_sha256']}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
