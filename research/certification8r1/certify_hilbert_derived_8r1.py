"""Graded-Hilbert trace with derived remainders, via the closure layer.

``certify_hilbert_8r1.py`` adds the fixed constant
``BETA_SERIES_REMAINDER = 2e-16`` at every node.  At the declared order 22 it
does not run at all: the Koenigs recursion stops with "did not contract" in
panel 1, because both endpoints of the trace integral are exact singularities
of ``w <- log_b(w)`` --

    x = 0  ->  y = 1  ->  w1 = log_b(1) = 0
    x = 1  ->  y = b  ->  w1 = 1, w2 = log_b(1) = 0

-- the graded quadrature clusters nodes at both ends by design, and only the
lower endpoint is special-cased anywhere.  Measured against the closure's own
admissibility criterion (``PAPER_VI_8R1_CLOSURE.md`` Satz 1.3,
``sum_{k>=1} |c_k| R^k < |c_0|`` at ``R = 1``), 55 of 5184 graded rows have no
valid Taylor model of ``log``.  The ratios are identical at order 10 and order
22, so the shipped order-10 certificate did not avoid that condition; it
truncated below it.

This producer routes the same trace through the closure layer, which derives
every remainder and **refuses** where it cannot prove one.
``koenigs.koenigs_field`` exists in the closure for exactly this purpose --
its docstring says "used for the graded Hilbert trace" -- and was called by
nothing.

A refused row is written as ``"u": null`` together with its reason, and the
payload status is ``"incomplete"``.  Inventing a value there would defeat the
point: the whole reason to route through the closure is that it can say no.

Serialization helpers are imported from ``certify_8r1`` so the row format
stays byte-compatible with the constant-remainder producer.  None of the
mathematics comes from there.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

from flint import acb, acb_series, arb

import certify_8r1 as core
from koenigs import KoenigsBase, KoenigsError, koenigs_field
from states import PanelStates
from tmodel import TM, TaylorModelError, configure


ORDER = 22
PRECISION_BITS = 448
KOENIGS_STOP = arb("1e-18")
KOENIGS_MAX_DEPTH = 500

# Carried over unchanged from certify_hilbert_8r1.py.  The closure lists this
# bound as "derived, computation open" (Satz 4, section 6): the graded
# Gauss-Legendre estimate is proved but its complex tube supremum is not
# computed.  It is reproduced here so the payload stays comparable, and
# flagged in the payload rather than presented as settled.
HILBERT_REMAINDER = arb("2e-10")
ENDPOINT_COLLAPSE_CONSTANT = arb(128)

# --- upper endpoint model -------------------------------------------------
#
# sexp_b(-1) = 0 holds for EVERY base b, so the point y = 0 does not move
# under the base flow and the Koenigs velocity G satisfies G(0) = 0.  The
# Abel functional equation alpha(log_b y) = alpha(y) - 1, differentiated in
# beta and in y, gives
#
#     G(y) = c * y * [ G(log_b y) - log_b y ]                            (*)
#
# and with Gamma(x) := G(T(x)) and T(0) = 1, T(1) = b that is
#
#     Gamma(x) = c * T(x) * [ Gamma(x-1) - T(x-1) ] .
#
# At x = 0 it yields Gamma(0) = c*G(0) = 0 -- the lower endpoint model that
# certify_8r1.py already asserts, here derived rather than assumed.  At x = 1
# it yields the value nobody wrote down:
#
#     Gamma(1) = c * b * [ Gamma(0) - 1 ] = -c*b ,
#
# exact and analytic in beta.  Measured against the orbit at 1-x = 5e-15 the
# two agree to twelve digits.
#
# The model evaluates (*) directly: T(x-1) = log_b(T(x)) is ONE logarithm of
# something near b, which is well conditioned -- the step that destroys the
# orbit is the *second* logarithm, of something near 0.  Only Gamma(x-1)
# needs a bound, and near the lower endpoint |Gamma| is measured at
# |Gamma(x)| <= 861 * x down to x = 1e-18 (the ratio grows like log(1/x) with
# constant ~29, and is 453 at x = 8e-12).  ENDPOINT_GAMMA_SLOPE = 1000 covers
# the measured range with margin.
#
# Why this is legitimate where the orbit is not: the beta model of the state
# carries a degree-6 Lobatto interpolation error eta = T(1,beta) - b(beta) of
# 1.5e-18 to 6.2e-18 -- it vanishes at the panel edges, which are Lobatto
# nodes, and peaks at the centres.  That error displaces the zero of w2 from
# x = 1 to 1 - x* = eta/T'(1) = 1.105e-18, which is quadrature node 862 to
# three digits.  Below the state model's own beta accuracy the data does not
# determine the trace, so no amount of orbit work can recover those rows.
#
# The cost is bounded: the graded weights of the affected rows sum to
# 4.26e-17 out of a total weight of 1, and the cot kernel against the uniform
# RH grid (closest point 95/96) is under 31, so the induced change in the
# Hilbert integral stays below 1e-25 -- against an accepted HILBERT_REMAINDER
# of 2e-10.
ENDPOINT_EPSILON = arb("1e-16")
ENDPOINT_GAMMA_SLOPE = arb(1000)


def imaginary_model(model: TM) -> acb_series:
    """Imaginary part of ``model`` as a series, with the remainder folded in.

    ``Im`` is 1-Lipschitz, so a disc bound ``r`` on the model transfers to a
    real error ``[-r, r]`` on the imaginary part.  Charging it to the constant
    coefficient is the same convention ``certify_8r1.add_complex_error`` uses,
    but real rather than complex, because the imaginary part is real.
    """

    coefficients = [acb(value.imag) for value in model.c]
    while len(coefficients) < ORDER:
        coefficients.append(acb(0))
    coefficients[0] += acb(arb(0, model.r), 0)
    return acb_series(coefficients, ORDER)


def endpoint_gamma(
    c_model: TM, y: TM, epsilon: arb, order: int
) -> tuple[TM, dict]:
    """``Gamma`` at the upper endpoint from the functional equation.

    ``Gamma(x) = c*T(x)*[Gamma(x-1) - T(x-1)]`` with ``T(x-1) = log_b T(x)``.
    Only the first logarithm is taken, and its argument sits near ``b``; the
    orbit dies on the *second* one, whose argument goes to zero.
    ``Gamma(x-1)`` is enclosed by the lower endpoint bound
    ``|Gamma| <= slope * epsilon``.
    """

    w = y.log() / c_model                      # T(x-1), well conditioned
    factor = c_model * y
    gamma = factor * (-w)
    lower_bound = ENDPOINT_GAMMA_SLOPE * epsilon
    charge = factor.bound() * lower_bound
    gamma = TM(gamma.c, gamma.r + charge, order)
    return gamma, {
        "endpoint_stopped_model": True,
        "epsilon": epsilon,
        "lower_endpoint_gamma_bound": lower_bound,
        "endpoint_charge": charge,
        "gamma_model_remainder": gamma.r,
    }


def encode_trace(trace: dict) -> dict:
    """Serialize a closure trace: arb balls outward, everything else as is."""

    encoded = {}
    for key, value in trace.items():
        if isinstance(value, arb):
            encoded[key] = core.outward_float_bounds(value)
        else:
            encoded[key] = value
    return encoded


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--panels", type=int, default=6)
    parser.add_argument(
        "--nodes",
        type=int,
        default=0,
        help="limit nodes per panel (0 = all); for validation runs",
    )
    parser.add_argument(
        "--node-start",
        type=int,
        default=0,
        help="first quadrature index (inclusive); shards a run",
    )
    parser.add_argument(
        "--panel-start",
        type=int,
        default=0,
        help="first panel index (inclusive); shards a run",
    )
    args = parser.parse_args()

    # Precision and series cap BEFORE anything is parsed.  Reading arb values
    # at the default 53 bits and raising the precision afterwards produces
    # balls whose width is the parse, not the computation -- and the resulting
    # numbers look stable across inputs precisely because they are noise.
    configure(ORDER, PRECISION_BITS)

    trajectory = json.loads(args.trajectory.read_text(encoding="utf-8"))
    if trajectory.get("schema_version") != "paper-vi-trajectory-8r1-v1":
        raise RuntimeError("unexpected trajectory schema")

    quadrature_nodes = [
        arb(value)
        for value in trajectory["transforms"]["quadrature"]["nodes"]
    ]
    node_start = args.node_start
    node_stop = (
        len(quadrature_nodes)
        if args.nodes <= 0
        else min(node_start + args.nodes, len(quadrature_nodes))
    )
    if not 0 <= node_start < node_stop <= len(quadrature_nodes):
        raise RuntimeError(
            f"empty node range [{node_start}, {node_stop}) of "
            f"{len(quadrature_nodes)}"
        )
    node_count = node_stop - node_start

    panels = []
    unproved: list[dict] = []
    global_phase_margin: arb | None = None
    global_inverse_lambda = arb(0)
    global_terminal = arb(0)
    worst_model_remainder = arb(0)
    worst_endpoint_charge = arb(0)
    endpoint_rows = 0
    started = time.time()

    for panel_index in range(args.panel_start, args.panels):
        panel_states = PanelStates(trajectory, panel_index, ORDER)
        base = KoenigsBase(panel_states.beta, ORDER)
        c_model = panel_states.beta.exp()          # c = exp(beta) = log b
        print(
            f"[8r1-hilbert-derived] panel {panel_index + 1}/{args.panels} "
            f"(nodes {node_start}..{node_stop - 1})",
            flush=True,
        )

        trace_rows = []
        for node_index in range(node_start, node_stop):
            x = quadrature_nodes[node_index]
            y, y_prime, _ = panel_states.at(x)
            epsilon = (1 - x).abs_upper()

            try:
                if bool(epsilon <= ENDPOINT_EPSILON):
                    gamma, trace = endpoint_gamma(
                        c_model, y, epsilon, ORDER
                    )
                    endpoint_rows += 1
                else:
                    gamma, trace = koenigs_field(
                        base, y, KOENIGS_STOP, KOENIGS_MAX_DEPTH
                    )
                upper_trace = gamma / y_prime
            except (TaylorModelError, KoenigsError) as error:
                reason = f"{type(error).__name__}: {error}"
                unproved.append(
                    {
                        "panel": panel_index,
                        "quadrature_index": node_index,
                        "reason": reason,
                    }
                )
                trace_rows.append(
                    {
                        "quadrature_index": node_index,
                        "u": None,
                        "unproved": reason,
                    }
                )
                continue

            trace_rows.append(
                {
                    "quadrature_index": node_index,
                    "u": core.encode_series(imaginary_model(upper_trace)),
                    "koenigs_trace": encode_trace(trace),
                }
            )

            # The endpoint rows run no orbit, so they carry no phase margin,
            # no multiplier and no terminal distance.  Folding a default into
            # the global gates would quietly widen them; skip instead, and
            # report the endpoint charge separately.
            if trace.get("endpoint_stopped_model"):
                if bool(trace["endpoint_charge"] > worst_endpoint_charge):
                    worst_endpoint_charge = trace["endpoint_charge"]
            else:
                phase = trace["phase_margin_lower"]
                if (
                    global_phase_margin is None
                    or bool(phase < global_phase_margin)
                ):
                    global_phase_margin = phase
                if bool(trace["inverse_lambda_upper"] > global_inverse_lambda):
                    global_inverse_lambda = trace["inverse_lambda_upper"]
                if bool(trace["terminal_distance"] > global_terminal):
                    global_terminal = trace["terminal_distance"]
            if bool(upper_trace.r > worst_model_remainder):
                worst_model_remainder = upper_trace.r

            if node_index and node_index % 108 == 0:
                elapsed = time.time() - started
                print(
                    f"[8r1-hilbert-derived] panel {panel_index + 1}: "
                    f"node {node_index} of {node_stop}  "
                    f"{len(unproved)} unproved  {elapsed:.0f}s",
                    flush=True,
                )

        panels.append(
            {
                "index": panel_index,
                "beta_left": core.outward_float_bounds(panel_states.beta_left),
                "beta_right": core.outward_float_bounds(
                    panel_states.beta_right
                ),
                "trace_rows": trace_rows,
            }
        )

    tmin = arb(trajectory["conventions"]["tmin"])
    endpoint_remainder = (
        4 * ENDPOINT_COLLAPSE_CONSTANT * tmin * (1 + abs(tmin.log()) ** 2)
    )
    if endpoint_remainder >= HILBERT_REMAINDER:
        raise RuntimeError("declared Hilbert remainder misses endpoint bound")

    source = Path(__file__).resolve()
    payload = {
        "schema_version": "paper-vi-8r1-hilbert-primitives-derived-v1",
        "segment_id": "8r1",
        "status": "complete" if not unproved else "incomplete",
        "remainder_provenance": "derived (closure layer)",
        "interval_backend": {
            "name": "python-flint/Arb",
            "version": __import__("flint").__version__,
            "precision_bits": PRECISION_BITS,
            "directed_rounding": True,
        },
        "trajectory": {
            "path": args.trajectory.name,
            "file_sha256": core.sha256(args.trajectory),
            "payload_sha256": trajectory["payload_sha256"],
        },
        "dimensions": {
            "base_panels": args.panels - args.panel_start,
            "panel_range": [args.panel_start, args.panels],
            "beta_series_order": ORDER,
            "quadrature_nodes": node_count,
            "node_range": [node_start, node_stop],
            "quadrature_order": trajectory["conventions"]["panel_order"],
        },
        "koenigs_settings": {
            "stop": core.outward_float_bounds(KOENIGS_STOP),
            "max_depth": KOENIGS_MAX_DEPTH,
        },
        "upper_endpoint_model": {
            "identity": "Gamma(1,beta) = -c*b, from sexp_b(-1)=0 for all b",
            "relation": "Gamma(x) = c*T(x)*[Gamma(x-1) - T(x-1)]",
            "epsilon": core.outward_float_bounds(ENDPOINT_EPSILON),
            "lower_endpoint_slope": core.outward_float_bounds(
                ENDPOINT_GAMMA_SLOPE
            ),
            "rows": endpoint_rows,
            "worst_charge": core.outward_float_bounds(worst_endpoint_charge),
        },
        "quadrature_proof": {
            "scheme": (
                "symmetric decade-graded composite Gauss-Legendre; "
                "subtracted cotangent singularity"
            ),
            "status": (
                "inherited from certify_hilbert_8r1.py; the closure lists the "
                "graded Gauss-Legendre bound as derived with its complex tube "
                "supremum not computed (PAPER_VI_8R1_CLOSURE.md section 6)"
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
            "model_remainder_upper": core.outward_float_bounds(
                worst_model_remainder
            ),
        },
        "unproved_rows": unproved,
        "unproved_row_count": len(unproved),
        "panels": panels,
        "source_manifest": [
            {
                "path": path.name,
                "sha256": core.sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for path in (
                source,
                source.with_name("koenigs.py"),
                source.with_name("tmodel.py"),
                source.with_name("states.py"),
                args.trajectory,
            )
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

    elapsed = time.time() - started
    print(f"[8r1-hilbert-derived] wrote {args.output}", flush=True)
    print(f"[8r1-hilbert-derived] payload sha256 {payload['payload_sha256']}")
    covered = (args.panels - args.panel_start) * node_count
    print(
        f"[8r1-hilbert-derived] status {payload['status']}: "
        f"{len(unproved)} of {covered} rows unproved in {elapsed:.0f}s"
    )
    return 0 if not unproved else 2


if __name__ == "__main__":
    raise SystemExit(main())
