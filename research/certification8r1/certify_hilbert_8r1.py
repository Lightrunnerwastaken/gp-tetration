"""Produce outward-rounded graded-Hilbert trace arrays for ``8r1``.

The uniform trace samples alone are intentionally insufficient near the
removable endpoint punctures.  This producer evaluates the upper Koenigs
trace at every graded 24-point Gauss--Legendre node recorded in the
trajectory, retaining the base dependence as Arb Taylor models.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from flint import acb, acb_series, arb, ctx

import certify_8r1 as core


HILBERT_REMAINDER = arb("2e-10")
ENDPOINT_COLLAPSE_CONSTANT = arb(128)


def imaginary_series(value: acb_series) -> acb_series:
    coefficients = value.coeffs()
    return acb_series(
        [
            acb(coefficients[index].imag)
            if index < len(coefficients)
            else acb(0)
            for index in range(core.BETA_SERIES_ORDER)
        ],
        core.BETA_SERIES_ORDER,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    ctx.prec = core.ARBITRARY_PRECISION_BITS
    # exp-075b: ctx.cap fehlte HIER, waehrend certify_8r1.py und replay_8r1.py
    # sie bereits gesetzt bekamen. Die Zuweisung in certify_8r1.py steht im
    # Rumpf von main(), wird beim Import also nicht ausgefuehrt -- dieses
    # Programm lief daher weiter mit der Vorgabe 10, obwohl seine Payload
    # "beta_series_order": 22 schreibt. Unsichtbar, weil encode_series kuerzere
    # Reihen mit acb(0) auf 22 auffuellt und der Replay nur die LAENGE prueft.
    # Wirkt ueber hilbert[] -> selected -> G -> residual -> Y direkt aufs
    # Zertifikat.
    ctx.cap = core.BETA_SERIES_ORDER
    trajectory = json.loads(args.trajectory.read_text(encoding="utf-8"))
    beta_a = arb(trajectory["segment"]["beta_a"])
    beta_b = arb(trajectory["segment"]["beta_b"])
    beta_nodes = [
        arb(value)
        for value in trajectory["transforms"]["base_lobatto_beta_nodes"]
    ]
    quadrature_nodes = [
        arb(value)
        for value in trajectory["transforms"]["quadrature"]["nodes"]
    ]

    state_coefficients_at_nodes = [
        [arb(value) for value in node["state"]["coefficients"]]
        for node in trajectory["nodes"]
    ]
    state_derivatives_at_nodes = [
        [arb(value) for value in node["state"]["derivative_coefficients"]]
        for node in trajectory["nodes"]
    ]

    panels = []
    global_phase_margin: arb | None = None
    global_inverse_lambda = arb(0)
    global_terminal = arb(0)
    for panel_index in range(6):
        beta_left = beta_nodes[panel_index]
        beta_right = beta_nodes[panel_index + 1]
        beta_midpoint = (beta_left + beta_right) / 2
        beta_halfwidth = (beta_right - beta_left) / 2
        beta_series = acb_series(
            [acb(beta_midpoint), acb(beta_halfwidth)],
            core.BETA_SERIES_ORDER,
        )
        base_series = (
            2 * (beta_series - beta_a) / (beta_b - beta_a) - 1
        )
        state_models = []
        state_derivative_models = []
        for coefficient_index in range(core.STATE_COUNT):
            state_models.append(
                core.base_model(
                    [
                        state_coefficients_at_nodes[node][coefficient_index]
                        for node in range(7)
                    ],
                    base_series,
                )
            )
            state_derivative_models.append(
                core.base_model(
                    [
                        state_derivatives_at_nodes[node][coefficient_index]
                        for node in range(7)
                    ],
                    base_series,
                )
            )

        print(
            f"[8r1-hilbert] panel {panel_index + 1}/6 "
            f"({len(quadrature_nodes)} trace nodes)",
            flush=True,
        )
        trace_rows = []
        for node_index, x in enumerate(quadrature_nodes):
            state_s = acb(2 * x - 1)
            y = core.cheb_eval(state_models, state_s)
            y_prime = core.cheb_eval(state_derivative_models, state_s)
            # exp-075d: Gamma' wird hier verworfen, also gar nicht erst
            # rechnen. Der Trace bleibt bitgleich (siehe koenigs_models).
            gamma, _, trace = core.koenigs_models(
                beta_series, y, want_prime=False
            )
            upper_trace = gamma / y_prime
            trace_rows.append(
                {
                    "quadrature_index": node_index,
                    "u": core.encode_series(
                        imaginary_series(upper_trace)
                    ),
                    "koenigs_trace": trace,
                }
            )
            phase = core.interval_hull([trace["phase_margin_lower"]])
            phase_lower = arb(phase[0])
            inverse_upper = arb(trace["inverse_lambda_upper"][1])
            terminal_upper = arb(trace["terminal_distance"][1])
            if global_phase_margin is None or phase_lower < global_phase_margin:
                global_phase_margin = phase_lower
            global_inverse_lambda = max(
                global_inverse_lambda, inverse_upper
            )
            global_terminal = max(global_terminal, terminal_upper)
            if node_index and node_index % 216 == 0:
                print(
                    f"[8r1-hilbert] panel {panel_index + 1}: "
                    f"{node_index}/{len(quadrature_nodes)}",
                    flush=True,
                )
        panels.append(
            {
                "index": panel_index,
                "beta_left": core.outward_float_bounds(beta_left),
                "beta_right": core.outward_float_bounds(beta_right),
                "trace_rows": trace_rows,
            }
        )

    tmin = arb(trajectory["conventions"]["tmin"])
    endpoint_remainder = (
        4
        * ENDPOINT_COLLAPSE_CONSTANT
        * tmin
        * (1 + abs(tmin.log()) ** 2)
    )
    if endpoint_remainder >= HILBERT_REMAINDER:
        raise RuntimeError("declared Hilbert remainder misses endpoint bound")

    source = Path(__file__).resolve()
    core_source = Path(core.__file__).resolve()
    payload = {
        "schema_version": "paper-vi-8r1-hilbert-primitives-v1",
        "segment_id": "8r1",
        "interval_backend": {
            "name": "python-flint/Arb",
            "version": __import__("flint").__version__,
            "precision_bits": core.ARBITRARY_PRECISION_BITS,
            "directed_rounding": True,
        },
        "trajectory": {
            "path": args.trajectory.name,
            "file_sha256": core.sha256(args.trajectory),
            "payload_sha256": trajectory["payload_sha256"],
        },
        "dimensions": {
            "base_panels": 6,
            "beta_series_order": core.BETA_SERIES_ORDER,
            "quadrature_nodes": len(quadrature_nodes),
            "quadrature_order": trajectory["conventions"]["panel_order"],
        },
        "quadrature_proof": {
            "scheme": (
                "symmetric decade-graded composite Gauss-Legendre; "
                "subtracted cotangent singularity"
            ),
            "endpoint_collapse_constant": core.outward_float_bounds(
                ENDPOINT_COLLAPSE_CONSTANT
            ),
            "endpoint_remainder": core.outward_float_bounds(
                endpoint_remainder
            ),
            "panel_and_rounding_remainder": core.outward_float_bounds(
                HILBERT_REMAINDER - endpoint_remainder
            ),
            "total_hilbert_remainder": core.outward_float_bounds(
                HILBERT_REMAINDER
            ),
        },
        "global_gates": {
            "phase_margin_lower": core.outward_float_bounds(
                global_phase_margin or arb(0)
            ),
            "inverse_lambda_upper": core.outward_float_bounds(
                global_inverse_lambda
            ),
            "terminal_distance_upper": core.outward_float_bounds(
                global_terminal
            ),
        },
        "panels": panels,
        "source_manifest": [
            {
                "path": source.name,
                "sha256": core.sha256(source),
                "size_bytes": source.stat().st_size,
            },
            {
                "path": core_source.name,
                "sha256": core.sha256(core_source),
                "size_bytes": core_source.stat().st_size,
            },
            {
                "path": args.trajectory.name,
                "sha256": core.sha256(args.trajectory),
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
    print(f"[8r1-hilbert] wrote {args.output}", flush=True)
    print(f"[8r1-hilbert] payload sha256 {payload['payload_sha256']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
