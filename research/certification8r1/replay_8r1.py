"""Independent array replay for the outward-rounded ``8r1`` primitives.

The replay deliberately does not import either numerical producer.  It reads
only serialized interval coefficients, reconstructs the RH/Hardy selection,
the complex tail, the normalized Volterra field, its residual and the
scale-Lipschitz majorant, and then chooses a uniform mathematical subdivision
of the physical 8r1 segment.  Every accepted scalar is recomputed with Arb.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from flint import acb, acb_poly, acb_series, arb, ctx


PRECISION_BITS = 320
# Two interval producers exist, with the same row layout.  The original
# charges every Koenigs and tail row the fixed constant BETA_SERIES_REMAINDER;
# certify_derived_8r1.py takes each remainder from the closure layer and
# writes a row as "unproved" where it cannot prove one.
PRIMITIVE_SCHEMAS = (
    "paper-vi-8r1-interval-primitives-v1",
    "paper-vi-8r1-interval-primitives-derived-v1",
)
RHO_STATE = arb("0.10")
RHO_ZERO = arb("0.05")
RHO_TRACE = arb("0.05")
RHO_WORK = arb("0.01")
RHO_END = arb("0.01")
SIGMA_GEOMETRIC = arb("0.10")
RESIDUAL_HEAD = 47
HUB_CUTOFF = 16
COMPLEX_DERIVATIVE_INFLATION = arb(16)
TARGET_Q_UPPER = arb("0.24")
RADIUS_SAFETY = arb("1.25")
M1_UPPER = arb("0.01")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def interval(pair: list[str]) -> arb:
    lower = arb(pair[0])
    upper = arb(pair[1])
    if lower > upper:
        raise RuntimeError(f"reversed interval endpoints: {pair}")
    midpoint = (lower + upper) / 2
    radius = (upper - lower) / 2
    return arb(midpoint, radius)


def complex_interval(record: dict) -> acb:
    return acb(interval(record["re"]), interval(record["im"]))


def series(record: list[dict], order: int) -> acb_series:
    if len(record) != order:
        raise RuntimeError("series coefficient count mismatch")
    return acb_series([complex_interval(value) for value in record], order)


def real_series(value: acb_series, order: int) -> acb_series:
    coefficients = value.coeffs()
    return acb_series(
        [
            acb(coefficients[index].real)
            if index < len(coefficients)
            else acb(0)
            for index in range(order)
        ],
        order,
    )


def imag_series(value: acb_series, order: int) -> acb_series:
    coefficients = value.coeffs()
    return acb_series(
        [
            acb(coefficients[index].imag)
            if index < len(coefficients)
            else acb(0)
            for index in range(order)
        ],
        order,
    )


def add_real_error(
    value: acb_series, error: arb, order: int
) -> acb_series:
    coefficients = value.coeffs()
    while len(coefficients) < order:
        coefficients.append(acb(0))
    coefficients[0] += acb(arb(0, error))
    return acb_series(coefficients, order)


def series_box(value: acb_series) -> acb:
    return acb_poly(value.coeffs())(arb("0 +/- 1"))


def dft(values: list[acb_series], order: int) -> dict[int, acb_series]:
    count = len(values)
    output: dict[int, acb_series] = {}
    for k in range(-count // 2 + 1, count // 2):
        total = acb_series([], order)
        for j, value in enumerate(values):
            angle = -2 * arb.pi() * k * j / count
            root = acb(angle.cos(), angle.sin())
            total += value * root
        output[k] = total / count
    return output


def inverse_dft(
    modes: dict[int, acb_series], count: int, order: int
) -> list[acb_series]:
    output = []
    for j in range(count):
        value = acb_series([], order)
        for k, mode in modes.items():
            angle = 2 * arb.pi() * k * j / count
            value += mode * acb(angle.cos(), angle.sin())
        output.append(value)
    return output


def abs_upper(value: acb | arb) -> arb:
    return value.abs_upper() if isinstance(value, acb) else value.abs_upper()


def min_abs_lower(values: list[acb_series]) -> arb:
    result: arb | None = None
    for value in values:
        candidate = series_box(value).abs_lower()
        if result is None or candidate < result:
            result = candidate
    assert result is not None
    return result


def max_abs_upper(values: list[acb_series]) -> arb:
    result = arb(0)
    for value in values:
        candidate = series_box(value).abs_upper()
        if candidate > result:
            result = candidate
    return result


def wiener_norm(
    modes: dict[int, acb_series], rho: arb, cutoff: int | None = None
) -> arb:
    result = arb(0)
    for k, value in modes.items():
        if cutoff is None or abs(k) <= cutoff:
            result += (
                series_box(value).abs_upper()
                * (2 * arb.pi() * rho * abs(k)).exp()
            )
    return result


def cheb_point(coefficients: list[arb], x: arb) -> arb:
    s = 2 * x - 1
    b1 = arb(0)
    b2 = arb(0)
    for coefficient in coefficients[:0:-1]:
        b0 = 2 * s * b1 - b2 + coefficient
        b2, b1 = b1, b0
    return s * b1 - b2 + coefficients[0]


def positive_margin(value: arb) -> bool:
    return bool(value > 0)


def upper_less(value: arb, bound: arb) -> bool:
    return bool(value < bound)


def interval_text(value: arb, digits: int = 18) -> str:
    return value.str(digits, radius=True)


def complex_text(value: acb) -> dict[str, str]:
    return {
        "re": interval_text(value.real),
        "im": interval_text(value.imag),
    }


def canonical_payload_hash(payload: dict) -> str:
    copy = dict(payload)
    copy.pop("payload_sha256", None)
    canonical = json.dumps(
        copy, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primitives", type=Path, required=True)
    parser.add_argument("--hilbert", type=Path, required=True)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    ctx.prec = PRECISION_BITS

    # exp-075 (2026-07-31): python-flint kappt JEDE Serienoperation bei der globalen
    # Laenge ctx.cap, deren Vorgabe 10 ist. Keines der drei Programme hat sie je
    # gesetzt, waehrend die Payload "beta_series_order": 22 behauptete -- alle
    # beta-Modelle waren in Wahrheit Polynome vom Grad 9, die Koeffizienten 10..21
    # exakt null. Gemeldet im Closure-Dokument, auf dieser Maschine unabhaengig
    # nachgeprueft (flint 0.9.0, ctx.cap == 10). Ohne diese Zeile erzeugt der Lauf
    # ein Artefakt, dessen Etikett nicht stimmt.
    primitives = json.loads(args.primitives.read_text(encoding="utf-8"))
    # ctx.cap MUSS aus den Primitiven kommen, nicht aus einer Konstante: der
    # Replay soll genau die Ordnung nachrechnen, die die Payload behauptet.
    # Stimmt sie nicht mit dem ueberein, was der Produzent konnte, faellt das
    # hier auf statt unbemerkt zu bleiben.
    ctx.cap = int(primitives["dimensions"]["beta_series_order"])
    hilbert_primitives = json.loads(
        args.hilbert.read_text(encoding="utf-8")
    )
    trajectory = json.loads(args.trajectory.read_text(encoding="utf-8"))
    issues: list[str] = []

    primitive_schema = primitives.get("schema_version")
    if primitive_schema not in PRIMITIVE_SCHEMAS:
        issues.append("unexpected primitive schema")
    if canonical_payload_hash(primitives) != primitives.get("payload_sha256"):
        issues.append("primitive canonical payload hash mismatch")
    # Same rule as for the Hilbert leg below: a derived payload with rows the
    # closure refused must not reach the arithmetic.  Every one of the 96
    # rows per panel enters the DFT, so a missing row is an undefined term.
    unproved_primitive_rows = primitives.get("unproved_row_count", 0)
    if unproved_primitive_rows:
        issues.append(
            f"primitive payload carries {unproved_primitive_rows} unproved "
            "rows; the interval arrays are undefined and the replay stops"
        )
    # Two Hilbert producers exist.  The original charges every row the fixed
    # constant BETA_SERIES_REMAINDER; the derived one takes each remainder
    # from the closure layer and writes "u": null where it cannot prove one.
    # Both are accepted here, but a payload carrying unproved rows must not
    # reach the quadrature: the graded Hilbert sum needs every row, so a
    # missing one is not a small gap to widen, it is an undefined term.
    hilbert_schema = hilbert_primitives.get("schema_version")
    if hilbert_schema not in (
        "paper-vi-8r1-hilbert-primitives-v1",
        "paper-vi-8r1-hilbert-primitives-derived-v1",
    ):
        issues.append("unexpected Hilbert primitive schema")
    unproved_rows = hilbert_primitives.get("unproved_row_count", 0)
    if unproved_rows:
        issues.append(
            f"Hilbert payload carries {unproved_rows} unproved rows; "
            "the graded Hilbert term is undefined and the replay stops"
        )
    if hilbert_primitives.get("short_panels"):
        issues.append(
            "Hilbert payload has panels with an incomplete node range: "
            f"{hilbert_primitives['short_panels']}"
        )
    if (
        canonical_payload_hash(hilbert_primitives)
        != hilbert_primitives.get("payload_sha256")
    ):
        issues.append("Hilbert canonical payload hash mismatch")
    if sha256(args.trajectory) != primitives["trajectory"]["file_sha256"]:
        issues.append("trajectory file hash mismatch")
    if (
        sha256(args.trajectory)
        != hilbert_primitives["trajectory"]["file_sha256"]
    ):
        issues.append("Hilbert trajectory file hash mismatch")
    if (
        trajectory.get("payload_sha256")
        != primitives["trajectory"]["payload_sha256"]
    ):
        issues.append("trajectory payload hash mismatch")
    if (
        trajectory.get("conventions_sha256")
        != primitives.get("conventions_sha256")
    ):
        issues.append("conventions hash mismatch")

    # Stop before the quadrature when a Hilbert row is missing.  Proceeding
    # would either crash on "u": null or, worse, invite someone to fill the
    # gap with a zero -- which reads as "this row contributes nothing" when
    # what it means is "nobody could bound this row".
    blocking = [
        issue
        for issue in issues
        if "unproved rows" in issue or "incomplete node range" in issue
    ]
    if blocking:
        report = {
            "schema_version": "paper-vi-8r1-replay-report-v1",
            "segment_id": "8r1",
            "passed": False,
            "guarantee": "replay rejected before quadrature",
            "issues": issues,
            "primitive_schema": primitive_schema,
            "hilbert_schema": hilbert_schema,
            "unproved_primitive_row_count": unproved_primitive_rows,
            "unproved_primitive_rows": primitives.get("unproved_rows", []),
            "unproved_row_count": unproved_rows,
            "unproved_rows": hilbert_primitives.get("unproved_rows", []),
            "source_hashes": {
                "hilbert_file_sha256": sha256(args.hilbert),
                "trajectory_file_sha256": sha256(args.trajectory),
                "replay_source_sha256": sha256(Path(__file__).resolve()),
            },
        }
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(
            json.dumps(report, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print("8r1 array replay: FAIL")
        for issue in issues:
            print(f"issue: {issue}")
        print(f"report: {args.report}")
        return 1

    # Re-evaluate normalization, endpoint functional equations, derivative
    # compatibility and the raw U=x+q+h transform from the trajectory arrays.
    normalization_gates = []
    max_normalization_defect = arb(0)
    max_functional_defect = arb(0)
    max_derivative_defect = arb(0)
    max_raw_transform_defect = arb(0)
    anchor_log_base = arb(
        trajectory["nodes"][0]["raw_transform"]["anchor_log_base"]
    )
    state_x_nodes = [
        arb(value)
        for value in trajectory["transforms"]["state_chebyshev_nodes"]
    ]
    for node in trajectory["nodes"]:
        coefficients = [arb(value) for value in node["state"]["coefficients"]]
        derivative_coefficients = [
            arb(value) for value in node["state"]["derivative_coefficients"]
        ]
        T0 = cheb_point(coefficients, arb(0))
        T1 = cheb_point(coefficients, arb(1))
        T0_prime = cheb_point(derivative_coefficients, arb(0))
        T1_prime = cheb_point(derivative_coefficients, arb(1))
        base = arb(node["base"])
        log_base = base.log()
        normalization_defect = abs(T0 - 1)
        functional_defect = abs(T1 - (log_base * T0).exp())
        derivative_defect = abs(
            T1_prime - log_base * T1 * T0_prime
        )
        alpha = arb(node["raw_transform"]["alpha"])
        offset = arb(node["raw_transform"]["offset"])
        transform_alpha_defect = abs(
            alpha - anchor_log_base / log_base
        )
        transform_offset_defect = abs(
            offset - alpha.log() / log_base
        )
        q_value = arb(node["raw_transform"]["q"])
        raw_U = [
            arb(value)
            for value in node["raw_transform"]["U_values_at_state_nodes"]
        ]
        raw_h = [
            arb(value)
            for value in node["normalized_hub"]["h_values_at_state_nodes"]
        ]
        raw_transform_defect = max(
            abs(
                raw_U[index]
                - state_x_nodes[index]
                - q_value
                - raw_h[index]
            )
            for index in range(len(state_x_nodes))
        )
        raw_transform_defect = max(
            raw_transform_defect,
            transform_alpha_defect,
            transform_offset_defect,
        )
        max_normalization_defect = max(
            max_normalization_defect, normalization_defect
        )
        max_functional_defect = max(
            max_functional_defect, functional_defect
        )
        max_derivative_defect = max(
            max_derivative_defect, derivative_defect
        )
        max_raw_transform_defect = max(
            max_raw_transform_defect, raw_transform_defect
        )
        normalization_gates.append(
            {
                "node": node["index"],
                "T0_minus_1": interval_text(normalization_defect),
                "T1_minus_exp_logb_T0": interval_text(
                    functional_defect
                ),
                "derivative_boundary_defect": interval_text(
                    derivative_defect
                ),
                "raw_transform_defect": interval_text(
                    raw_transform_defect
                ),
            }
        )
    representation_gate = arb("1e-44")
    if max_normalization_defect >= representation_gate:
        issues.append("normalization representation gate failed")
    if max_functional_defect >= representation_gate:
        issues.append("functional-equation representation gate failed")
    if max_derivative_defect >= representation_gate:
        issues.append("derivative functional-equation gate failed")
    if max_raw_transform_defect >= representation_gate:
        issues.append("raw transformation gate failed")

    order = primitives["dimensions"]["beta_series_order"]
    count = primitives["dimensions"]["rh_uniform_count"]
    if count != 96 or len(primitives.get("panels", [])) != 6:
        issues.append("unexpected replay dimensions")
    if len(hilbert_primitives.get("panels", [])) != 6:
        issues.append("unexpected Hilbert replay dimensions")

    quadrature_nodes = [
        arb(value)
        for value in trajectory["transforms"]["quadrature"]["nodes"]
    ]
    quadrature_weights = [
        arb(value)
        for value in trajectory["transforms"]["quadrature"]["weights"]
    ]
    hilbert_remainder = interval(
        hilbert_primitives["quadrature_proof"]["total_hilbert_remainder"]
    ).upper()

    panel_results = []
    global_residual_head = arb(0)
    global_h_norm_state = arb(0)
    global_h_prime_state = arb(0)
    global_h_prime_trace = arb(0)
    global_m_t_prime: arb | None = None
    global_m_u_prime: arb | None = None
    global_m_q_prime: arb | None = None
    global_M_t = arb(0)
    global_M_t_prime = arb(0)
    global_M_t_second = arb(0)
    global_M_h_second = arb(0)
    global_M_q_prime = arb(0)
    global_M_q_second = arb(0)
    global_M_gamma = arb(0)
    global_M_gamma_prime = arb(0)
    global_M_tail = arb(0)
    global_M_tail_T = arb(0)
    global_M_tail_T_prime = arb(0)
    global_M_G_head = arb(0)
    global_T_real_lower: arb | None = None
    global_T_real_upper: arb | None = None
    global_tail_last_upper = arb(0)
    global_tail_ratio_upper = arb(0)
    global_tail_remainder_upper = arb(0)
    derived_primitives = (
        primitive_schema == "paper-vi-8r1-interval-primitives-derived-v1"
    )

    for panel, hilbert_panel in zip(
        primitives["panels"], hilbert_primitives["panels"]
    ):
        if panel["index"] != hilbert_panel["index"]:
            raise RuntimeError("base-panel mismatch between primitive files")
        rows = []
        for source in panel["primitive_rows"]:
            tail_last = interval(
                source["tail_trace"]["last_term_upper"]
            ).upper()
            global_tail_last_upper = max(
                global_tail_last_upper, tail_last
            )
            if derived_primitives:
                global_tail_ratio_upper = max(
                    global_tail_ratio_upper,
                    interval(source["tail_trace"]["ratio_upper"]).upper(),
                )
                global_tail_remainder_upper = max(
                    global_tail_remainder_upper,
                    interval(source["tail_trace"]["remainder"]).upper(),
                )
            rows.append(
                {
                    key: series(source[key], order)
                    for key in (
                        "T",
                        "T_prime",
                        "T_second",
                        "h",
                        "h_prime",
                        "h_second",
                        "h_beta",
                        "Gamma",
                        "Gamma_prime",
                        "tail",
                        "tail_partial_T",
                        "tail_partial_T_prime",
                    )
                }
            )

        T = [row["T"] for row in rows]
        T_prime = [row["T_prime"] for row in rows]
        T_second = [row["T_second"] for row in rows]
        h = [row["h"] for row in rows]
        h_prime = [row["h_prime"] for row in rows]
        h_second = [row["h_second"] for row in rows]
        h_beta = [row["h_beta"] for row in rows]
        Gamma = [row["Gamma"] for row in rows]
        Gamma_prime = [row["Gamma_prime"] for row in rows]
        tails = [row["tail"] for row in rows]
        tail_T = [row["tail_partial_T"] for row in rows]
        tail_T_prime = [row["tail_partial_T_prime"] for row in rows]

        v_geo = [Gamma[j] / T_prime[j] for j in range(count)]
        u = [imag_series(value, order) for value in v_geo]
        real_v = [real_series(value, order) for value in v_geo]
        quadrature_u = [
            series(source["u"], order)
            for source in hilbert_panel["trace_rows"]
        ]
        if len(quadrature_u) != len(quadrature_nodes):
            raise RuntimeError("graded Hilbert row count mismatch")
        hilbert = []
        for x_index in range(count):
            x = arb(x_index) / count
            value = acb_series([], order)
            for t, weight, u_t in zip(
                quadrature_nodes, quadrature_weights, quadrature_u
            ):
                kernel = (arb.pi() * (x - t)).cot()
                value += weight * (u_t - u[x_index]) * kernel
            hilbert.append(
                add_real_error(value, hilbert_remainder, order)
            )
        selected = [real_v[j] + hilbert[j] for j in range(count)]
        selected_zero = selected[0]
        selected = [value - selected_zero for value in selected]
        G = [selected[j] - tails[j] for j in range(count)]
        raw_rhs = [-(1 + h_prime[j]) * G[j] for j in range(count)]
        rhs_zero = raw_rhs[0]
        rhs = [value - rhs_zero for value in raw_rhs]
        residual = [h_beta[j] - rhs[j] for j in range(count)]

        h_modes = dft(h, order)
        h_prime_modes = dft(h_prime, order)
        selected_modes = dft(selected, order)
        G_modes = dft(G, order)
        residual_modes = dft(residual, order)
        residual_head = wiener_norm(
            residual_modes, RHO_WORK, RESIDUAL_HEAD
        )
        h_norm_state = wiener_norm(h_modes, RHO_STATE)
        h_prime_state = wiener_norm(h_prime_modes, RHO_STATE)
        h_prime_trace = wiener_norm(h_prime_modes, RHO_TRACE)
        G_head = wiener_norm(G_modes, RHO_WORK, RESIDUAL_HEAD)

        m_t_prime = min_abs_lower(T_prime)
        M_t = max_abs_upper(T)
        M_t_prime = max_abs_upper(T_prime)
        M_t_second = max_abs_upper(T_second)
        M_h_second = max_abs_upper(h_second)
        M_gamma = max_abs_upper(Gamma)
        M_gamma_prime = max_abs_upper(Gamma_prime)
        M_tail = max_abs_upper(tails)
        M_tail_T = max_abs_upper(tail_T)
        M_tail_T_prime = max_abs_upper(tail_T_prime)
        for value in T:
            value_box = series_box(value)
            real_lower = value_box.real.lower()
            real_upper = value_box.real.upper()
            global_T_real_lower = (
                real_lower
                if global_T_real_lower is None
                else min(global_T_real_lower, real_lower)
            )
            global_T_real_upper = (
                real_upper
                if global_T_real_upper is None
                else max(global_T_real_upper, real_upper)
            )

        m_u_prime = arb(1) - h_prime_state
        q_prime_models = [
            T_prime[j] / (1 + h_prime[j]) for j in range(count)
        ]
        q_second_models = [
            (
                T_second[j] * (1 + h_prime[j])
                - T_prime[j] * h_second[j]
            )
            / (1 + h_prime[j]) ** 3
            for j in range(count)
        ]
        m_q_prime = min_abs_lower(q_prime_models)
        M_q_prime = max_abs_upper(q_prime_models)
        M_q_second = max_abs_upper(q_second_models)

        global_residual_head = max(global_residual_head, residual_head)
        global_h_norm_state = max(global_h_norm_state, h_norm_state)
        global_h_prime_state = max(
            global_h_prime_state, h_prime_state
        )
        global_h_prime_trace = max(
            global_h_prime_trace, h_prime_trace
        )
        global_m_t_prime = (
            m_t_prime
            if global_m_t_prime is None
            else min(global_m_t_prime, m_t_prime)
        )
        global_m_u_prime = (
            m_u_prime
            if global_m_u_prime is None
            else min(global_m_u_prime, m_u_prime)
        )
        global_m_q_prime = (
            m_q_prime
            if global_m_q_prime is None
            else min(global_m_q_prime, m_q_prime)
        )
        global_M_t = max(global_M_t, M_t)
        global_M_t_prime = max(global_M_t_prime, M_t_prime)
        global_M_t_second = max(global_M_t_second, M_t_second)
        global_M_h_second = max(global_M_h_second, M_h_second)
        global_M_q_prime = max(global_M_q_prime, M_q_prime)
        global_M_q_second = max(global_M_q_second, M_q_second)
        global_M_gamma = max(global_M_gamma, M_gamma)
        global_M_gamma_prime = max(
            global_M_gamma_prime, M_gamma_prime
        )
        global_M_tail = max(global_M_tail, M_tail)
        global_M_tail_T = max(global_M_tail_T, M_tail_T)
        global_M_tail_T_prime = max(
            global_M_tail_T_prime, M_tail_T_prime
        )
        global_M_G_head = max(global_M_G_head, G_head)

        panel_results.append(
            {
                "index": panel["index"],
                "residual_head": interval_text(residual_head),
                "h_norm_rho_state": interval_text(h_norm_state),
                "h_prime_norm_rho_state": interval_text(h_prime_state),
                "m_T_prime_samples": interval_text(m_t_prime),
                "m_U_prime": interval_text(m_u_prime),
                "m_Q_prime_samples": interval_text(m_q_prime),
                "M_G_head": interval_text(G_head),
                "selected_RH_modes": [
                    {
                        "k": k,
                        **complex_text(series_box(selected_modes[k])),
                    }
                    for k in range(-47, 48)
                ],
                "kernel_G_modes": [
                    {
                        "k": k,
                        **complex_text(series_box(G_modes[k])),
                    }
                    for k in range(-47, 48)
                ],
                "volterra_residual_modes": [
                    {
                        "k": k,
                        **complex_text(series_box(residual_modes[k])),
                    }
                    for k in range(-47, 48)
                ],
            }
        )

    assert global_m_t_prime is not None
    assert global_m_u_prime is not None
    assert global_m_q_prime is not None
    assert global_T_real_lower is not None
    assert global_T_real_upper is not None

    # Extend the real interval cover to the stored nonperiodic tube.  The
    # inflation is deliberately explicit and is checked against the phase,
    # cut and zero margins below.
    tube_width = SIGMA_GEOMETRIC - RHO_TRACE
    M_gamma_prime_tube = (
        COMPLEX_DERIVATIVE_INFLATION * global_M_gamma_prime
    )
    M_gamma_tube = (
        global_M_gamma
        + tube_width * M_gamma_prime_tube * global_M_t_prime
    )
    m_t_prime_tube = (
        global_m_t_prime - tube_width * global_M_t_second
    )
    M_tail_tube = (
        global_M_tail
        + tube_width
        * (
            global_M_tail_T * global_M_t_prime
            + global_M_tail_T_prime * global_M_t_second
        )
        * COMPLEX_DERIVATIVE_INFLATION
    )
    if m_t_prime_tube <= 0:
        issues.append("T' tube exclusion failed")

    c_per = (
        1 + (-2 * arb.pi() * (SIGMA_GEOMETRIC - RHO_TRACE)).exp()
    ) / (
        1 - (-2 * arb.pi() * (SIGMA_GEOMETRIC - RHO_TRACE)).exp()
    )
    M_v_tube = M_gamma_tube / m_t_prime_tube
    M_R_trace = (2 + 4 * c_per) * M_v_tube
    M_G_trace = M_R_trace + M_tail_tube

    trace_gap = RHO_TRACE - RHO_WORK
    tail_factor_47 = (
        (-2 * arb.pi() * trace_gap * (RESIDUAL_HEAD + 1)).exp()
    )
    tail_factor_31 = (
        (-2 * arb.pi() * trace_gap * 32).exp()
    )
    alias_factor = 2 * (
        -2
        * arb.pi()
        * trace_gap
        * (count - RESIDUAL_HEAD)
    ).exp()
    G_tail_47 = M_G_trace * (tail_factor_47 + alias_factor)
    G_tail_31 = M_G_trace * tail_factor_31
    residual_tail = (
        (1 + global_h_prime_trace) * G_tail_47
        + global_h_prime_trace * G_tail_31
    )
    residual_total = global_residual_head + residual_tail

    # Paper-VI scale-Lipschitz majorant, reconstructed from primitive arrays.
    C_T = RHO_STATE * global_M_q_prime
    C_T_prime = (
        global_M_q_prime
        + RHO_STATE
        * global_M_q_second
        * (1 + global_h_prime_state)
    )
    C_geo = (
        M_gamma_prime_tube * C_T / m_t_prime_tube
        + M_gamma_tube * C_T_prime / (m_t_prime_tube**2)
    )
    C_trace = c_per * C_geo
    C_RH = 2 * (C_geo + 2 * C_trace)
    C_tail = (
        global_M_tail_T * global_M_q_prime * RHO_STATE
        + global_M_tail_T_prime
        * (
            global_M_q_prime
            + RHO_STATE
            * global_M_q_second
            * (1 + global_h_prime_state)
        )
    )
    C_G = C_RH + C_tail
    L = 2 * (
        M_G_trace + (1 + global_h_prime_state) * C_G
    )

    beta_a = arb(trajectory["segment"]["beta_a"])
    beta_b = arb(trajectory["segment"]["beta_b"])
    total_length = abs(beta_b - beta_a)
    microsegments = 1
    while True:
        cell_length = total_length / microsegments
        kappa = (
            arb("0.75")
            * (RHO_ZERO - RHO_END)
            / cell_length
        )
        q = 4 * arb(2).sqrt() * L / kappa
        if q < TARGET_Q_UPPER:
            break
        microsegments *= 2
        if microsegments > 2**30:
            issues.append("subdivision limit exceeded")
            break

    strip_slack = RHO_ZERO - kappa * cell_length - RHO_END
    Y_head = cell_length * RHO_ZERO.sqrt() * global_residual_head
    Y_tail = cell_length * RHO_ZERO.sqrt() * residual_tail
    Y = Y_head + Y_tail
    radius = RADIUS_SAFETY * Y / (1 - q)
    polynomial = Y + (q - 1) * radius
    endpoint_per_cell = radius / strip_slack.sqrt()
    endpoint_total = (
        microsegments
        * endpoint_per_cell
        * (M1_UPPER * total_length).exp()
    )

    # Full-ball Z1--Z7 margins.
    z1_margin = arb("0.10") - global_h_norm_state - endpoint_total
    derivative_ball = endpoint_total / (RHO_STATE - RHO_ZERO)
    z2_u_margin = global_m_u_prime - derivative_ball
    z2_q_margin = global_m_q_prime - derivative_ball
    z3_t_margin = m_t_prime_tube - derivative_ball
    # exp-075c: genau im Fehlerfall stuerzte dieser Block ab. Bei Z3-Versagen
    # ist z3_t_margin <= 0, damit wird 1/z3_t_margin negativ oder unendlich und
    # math.log2 wirft (domain error bzw. OverflowError bei 2**ceil(inf)). Der
    # Report entsteht aber erst am Ende von main(), also gab es bei Z3-Versagen
    # KEINEN Report, nur einen Traceback -- und damit war kein Mutationstest
    # gegen den Replay moeglich: ein FAIL war von einem Absturz nicht zu
    # unterscheiden. Die Pruefung m_t_prime_tube <= 0 weiter oben meldet das
    # Problem zwar, der Code rechnete danach trotzdem weiter.
    # Bei Versagen werden beide Groessen auf null gesetzt: das Z3-Urteil weiter
    # unten lautet dann korrekt False, statt dass der Lauf verschwindet.
    if z3_t_margin <= 0:
        issues.append("Z3: T' tube margin non-positive; reciprocal bound skipped")
        reciprocal_bound = arb(0)
        z3_reciprocal_margin = arb(0)
    else:
        reciprocal_needed = 1 / z3_t_margin
        reciprocal_mid = float(reciprocal_needed.upper())
        if not math.isfinite(reciprocal_mid) or reciprocal_mid <= 0:
            issues.append("Z3: reciprocal bound not finite")
            reciprocal_bound = arb(0)
            z3_reciprocal_margin = arb(0)
        else:
            reciprocal_bound = arb(2 ** math.ceil(math.log2(reciprocal_mid)))
            z3_reciprocal_margin = reciprocal_bound - reciprocal_needed

    beta_mid = (beta_a + beta_b) / 2
    beta_rad = total_length / 2
    beta_box = arb(beta_mid, beta_rad)
    c_box = beta_box.exp()
    L_fixed = -acb(-c_box).lambertw(-1) / c_box
    lambda_box = c_box * L_fixed
    z6_one_minus = abs(1 - lambda_box).abs_lower()
    z6_expansion = lambda_box.abs_lower() - 1

    # Replay a tube enclosure for T_0,T_1,T_2.  This supplies both the
    # logarithmic-cut distance (positive real part) and the a_j lower bounds
    # used by the finite/tail split.  The superexponential continuation is
    # checked separately from the stored reciprocal-product tail.
    real_mid = (global_T_real_lower + global_T_real_upper) / 2
    real_rad = (global_T_real_upper - global_T_real_lower) / 2
    imaginary_rad = (
        tube_width * global_M_t_prime + endpoint_total
    )
    transfer = acb(arb(real_mid, real_rad), arb(0, imaginary_rad))
    cut_margins = []
    transfer_margins = []
    fixed_point_margins = []
    for _ in range(3):
        cut_margins.append(transfer.real.lower())
        transfer_margins.append(transfer.abs_lower())
        fixed_point_margins.append((transfer - L_fixed).abs_lower())
        transfer = (c_box * transfer).exp()
    z4_cut_margin = (
        min(cut_margins + fixed_point_margins) - endpoint_total
    )
    tail_stop_lower = interval(
        primitives["proof_constants"]["tail_stop"]
    ).lower()
    if derived_primitives:
        # The closure's tail does not run until the last term drops below
        # TAIL_STOP: it stops as soon as the geometric domination
        # |a_{j+1}/a_j| = 1/(c y_{j+1}) is certified, and charges the proved
        # remainder last_term * ratio/(1 - ratio) into the array itself
        # (tails.certified_tail).  Its convergence statement is therefore
        # "ratio < 1 on every row", not "last term < stop"; measuring the
        # last term against TAIL_STOP here would reject a sound payload for
        # having stopped earlier than the constant producer.
        tail_convergence_rule = (
            "certified ratio rule 1/(c y_{j+1}) < 1 with checked forward "
            "invariance; remainder folded into the tail arrays"
        )
        tail_convergence_margin = 1 - global_tail_ratio_upper
    else:
        tail_convergence_rule = "last term below tail_stop"
        tail_convergence_margin = (
            tail_stop_lower - global_tail_last_upper
        )
    z5_transfer_margin = min(
        min(transfer_margins) - endpoint_total,
        tail_convergence_margin,
    )

    phase_margin = interval(
        primitives["global_koenigs_gates"]["phase_margin_hull"]
    ).lower()
    inverse_lambda = interval(
        primitives["global_koenigs_gates"]["inverse_lambda_hull"]
    ).upper()
    koenigs_contraction_margin = 1 - inverse_lambda
    z7_margin = min(
        phase_margin,
        koenigs_contraction_margin,
        z3_t_margin,
        arb("0.25") - endpoint_total,
    )

    margins = {
        "Z1_hub_ball": z1_margin,
        "Z2_U_prime": z2_u_margin,
        "Z2_Q_prime": z2_q_margin,
        "Z3_T_prime": z3_t_margin,
        "Z3_reciprocal": z3_reciprocal_margin,
        "Z4_cut_distance": z4_cut_margin,
        "Z5_transfer_tail": z5_transfer_margin,
        "Z6_one_minus_lambda": z6_one_minus,
        "Z6_expansion": z6_expansion,
        "Z7_branch_overlap": z7_margin,
    }
    for name, margin in margins.items():
        if not positive_margin(margin):
            issues.append(f"{name} is not strictly positive")
    if not (q >= 0 and q < 1):
        issues.append("q is not in [0,1)")
    if not polynomial < 0:
        issues.append("radii polynomial is not strictly negative")
    if strip_slack <= 0:
        issues.append("endpoint strip slack is not positive")

    checks = {
        "hash_bindings": not any("hash" in issue for issue in issues),
        "q_lt_1": bool(q < 1),
        "p_lt_0": bool(polynomial < 0),
        "strip_slack_positive": bool(strip_slack > 0),
        "Z1": bool(z1_margin > 0),
        "Z2": bool(z2_u_margin > 0 and z2_q_margin > 0),
        "Z3": bool(z3_t_margin > 0 and z3_reciprocal_margin > 0),
        "Z4": bool(z4_cut_margin > 0),
        "Z5": bool(z5_transfer_margin > 0),
        "Z6": bool(z6_one_minus > 0 and z6_expansion > 0),
        "Z7": bool(z7_margin > 0),
        "normalization": bool(
            max_normalization_defect < representation_gate
        ),
        "functional_equation": bool(
            max_functional_defect < representation_gate
            and max_derivative_defect < representation_gate
        ),
        "raw_transform": bool(
            max_raw_transform_defect < representation_gate
        ),
    }
    passed = not issues and all(checks.values())

    report = {
        "schema_version": "paper-vi-8r1-array-replay-v1",
        "segment_id": "8r1",
        "passed": passed,
        "guarantee": (
            "outward-rounded composite Ovsyannikov segment certificate"
            if passed
            else "replay rejected"
        ),
        "interval_backend": {
            "name": "python-flint/Arb",
            "precision_bits": PRECISION_BITS,
            "directed_rounding": True,
        },
        # Which producer each leg came from, so a reader can tell a mixed
        # verdict from an all-derived one without opening the payloads.
        "provenance": {
            "primitive_schema": primitive_schema,
            "primitive_remainders": primitives.get(
                "remainder_provenance",
                "set constant (certify_8r1.py, BETA_SERIES_REMAINDER)",
            ),
            "hilbert_schema": hilbert_schema,
            "hilbert_remainders": hilbert_primitives.get(
                "remainder_provenance",
                "set constant (certify_hilbert_8r1.py, BETA_SERIES_REMAINDER)",
            ),
        },
        "source_hashes": {
            "primitives_file_sha256": sha256(args.primitives),
            "primitives_payload_sha256": primitives["payload_sha256"],
            "hilbert_file_sha256": sha256(args.hilbert),
            "hilbert_payload_sha256": hilbert_primitives["payload_sha256"],
            "trajectory_file_sha256": sha256(args.trajectory),
            "trajectory_payload_sha256": trajectory["payload_sha256"],
            "conventions_sha256": trajectory["conventions_sha256"],
            "replay_source_sha256": sha256(Path(__file__).resolve()),
        },
        "segment": {
            "beta_a": trajectory["segment"]["beta_a"],
            "beta_b": trajectory["segment"]["beta_b"],
            "base_a": trajectory["segment"]["base_a"],
            "base_b": trajectory["segment"]["base_b"],
            "microsegment_count": microsegments,
            "cell_length": interval_text(cell_length),
            "rho0": interval_text(RHO_ZERO),
            "rho_end": interval_text(RHO_END),
            "kappa": interval_text(kappa),
            "strip_slack": interval_text(strip_slack),
        },
        "radii_intervals": {
            "Y_head": interval_text(Y_head),
            "Y_tail": interval_text(Y_tail),
            "Y": interval_text(Y),
            "L": interval_text(L),
            "q": interval_text(q),
            "r": interval_text(radius),
            "p(r)": interval_text(polynomial),
            "endpoint_per_cell": interval_text(endpoint_per_cell),
            "endpoint_total": interval_text(endpoint_total),
        },
        "operator_bounds": {
            "residual_head": interval_text(global_residual_head),
            "residual_tail": interval_text(residual_tail),
            "residual_total": interval_text(residual_total),
            "m_T_prime_tube": interval_text(m_t_prime_tube),
            "mu_T": interval_text(reciprocal_bound),
            "M_G_trace": interval_text(M_G_trace),
            "C_per": interval_text(c_per),
            "C_geo": interval_text(C_geo),
            "C_trace": interval_text(C_trace),
            "C_RH": interval_text(C_RH),
            "C_tail": interval_text(C_tail),
            "C_G": interval_text(C_G),
        },
        "normalization_and_transform_replay": {
            "representation_gate": interval_text(representation_gate),
            "max_normalization_defect": interval_text(
                max_normalization_defect
            ),
            "max_functional_defect": interval_text(
                max_functional_defect
            ),
            "max_derivative_defect": interval_text(
                max_derivative_defect
            ),
            "max_raw_transform_defect": interval_text(
                max_raw_transform_defect
            ),
            "nodes": normalization_gates,
        },
        "margins": {
            name: interval_text(value) for name, value in margins.items()
        },
        "transfer_replay": {
            "J0": 2,
            "cut_margins": [interval_text(value) for value in cut_margins],
            "fixed_point_margins": [
                interval_text(value) for value in fixed_point_margins
            ],
            "a_j": [
                interval_text(value) for value in transfer_margins
            ],
            "tail_last_term": interval_text(global_tail_last_upper),
            "tail_stop": interval_text(tail_stop_lower),
            "tail_convergence_rule": tail_convergence_rule,
            "tail_ratio_upper": (
                interval_text(global_tail_ratio_upper)
                if derived_primitives
                else None
            ),
            "tail_remainder_upper": (
                interval_text(global_tail_remainder_upper)
                if derived_primitives
                else None
            ),
            "tail_convergence_margin": interval_text(
                tail_convergence_margin
            ),
        },
        "panel_replay": panel_results,
        "checks": checks,
        "issues": issues,
    }
    canonical = json.dumps(
        report, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    report["payload_sha256"] = hashlib.sha256(canonical).hexdigest()
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print(f"8r1 array replay: {'PASS' if passed else 'FAIL'}")
    for name in ("Y", "L", "q", "r", "p(r)"):
        print(f"{name} = {report['radii_intervals'][name]}")
    for name in ("Z1", "Z2", "Z3", "Z4", "Z5", "Z6", "Z7"):
        print(f"{name}: {'PASS' if checks[name] else 'FAIL'}")
    if issues:
        for issue in issues:
            print(f"issue: {issue}")
    print(f"report: {args.report}")
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
