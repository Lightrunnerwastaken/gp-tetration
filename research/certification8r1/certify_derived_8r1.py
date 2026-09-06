"""Interval primitives with *derived* remainders, via the closure layer.

``certify_8r1.py`` builds the Koenigs, quotient and tail arrays of segment
``8r1`` and charges every one of them the fixed constant
``BETA_SERIES_REMAINDER = 2e-16`` -- the same kind of set constant that the
closure layer (``tmodel.py``, ``koenigs.py``, ``tails.py``) was written to
eliminate, and that ``certify_hilbert_derived_8r1.py`` already removed from
the graded Hilbert leg.  As long as only one leg went through the closure,
the replay verdict was of *mixed provenance*.

This producer routes the other leg.  It writes the **same row layout** as the
constant-remainder producer -- ``T``, ``T_prime``, ``T_second``, ``h``,
``h_prime``, ``h_second``, ``h_beta``, ``Gamma``, ``Gamma_prime``, ``tail``,
``tail_partial_T``, ``tail_partial_T_prime`` -- so that ``replay_8r1.py``
consumes it unchanged, but every transcendental remainder is derived:

* ``Gamma`` and ``Gamma_prime`` come from :func:`koenigs.koenigs_pair`, whose
  stopping error is telescoped from quantities the run itself produces (the
  asserted factors 128 and 1024 are not used), and whose ``Gamma'`` error is
  a Cauchy estimate on a closed state disc;
* the tails come from :func:`tails.certified_tail`, whose truncation rule is
  proved (ratio identity plus a checked forward-invariance premise) instead of
  the assumed ``R <= 2|a_N|``.

The state arrays ``T``, ``T'``, ``T''`` and the hub arrays ``h``, ``h'``,
``h''``, ``h_beta`` involve no transcendental step: they are Chebyshev
interpolation in the base and an exact DFT of stored values, so they carry no
remainder in either producer and are built with the same helpers.

Where the closure cannot prove a remainder it **refuses**; the row is written
as ``"unproved"`` with the reason, the payload status becomes ``incomplete``,
and ``replay_8r1.py`` stops before using the arrays.  Inventing a value there
would defeat the purpose of routing through the closure.

Serialization helpers are imported from ``certify_8r1`` so the row format
stays byte-compatible.  None of the remainders come from there.

Run time is dominated by the two stopped orbits per node (thin and fattened)
that ``koenigs_pair`` needs; the six base panels are independent and run as
separate processes (``--jobs``).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from flint import acb, acb_series, arb

import certify_8r1 as core
from koenigs import KoenigsBase, KoenigsError, koenigs_pair
from states import PanelStates
from tails import TailError, certified_tail
from tmodel import TM, TaylorModelError, configure


ORDER = 22
PRECISION_BITS = 448
STATE_COUNT = 48
RH_COUNT = 96
PANEL_COUNT = 6
KOENIGS_MAX_DEPTH = 500
# Thresholds are kept as text and parsed only after configure(): an arb
# parsed at the import-time default of 53 bits is a ball whose width is the
# parse, not the mathematics.
KOENIGS_STOP_TEXT = "1e-18"
STATE_DISC_RADIUS_TEXT = "1e-3"
TAIL_STOP_TEXT = "1e-50"
SCHEMA = "paper-vi-8r1-interval-primitives-derived-v1"


def encode_model(model: TM) -> list[dict[str, list[str]]]:
    """Serialize a Taylor model in the constant producer's row format.

    A model ``(c, r)`` encloses ``f(s)`` in ``P(s) + D(0, r)`` for every
    ``|s| <= 1``.  Folding ``r`` into the constant coefficient as a complex
    ball ``[-r, r] + i[-r, r]`` is the convention of
    ``certify_8r1.add_complex_error``; the square contains the disc, so the
    serialized row encloses the model.
    """

    coefficients = list(model.c)
    coefficients[0] += acb(arb(0, model.r), arb(0, model.r))
    return core.encode_series(acb_series(coefficients, ORDER))


def encode_trace(trace: dict) -> dict:
    return {
        key: core.outward_float_bounds(value) if isinstance(value, arb) else value
        for key, value in trace.items()
    }


def produce_panel(task: tuple[str, int, int]) -> dict:
    """Build one base panel.  Runs in its own process."""

    trajectory_path, panel_index, node_count = task
    configure(ORDER, PRECISION_BITS)
    koenigs_stop = arb(KOENIGS_STOP_TEXT)
    state_radius = arb(STATE_DISC_RADIUS_TEXT)
    tail_stop = arb(TAIL_STOP_TEXT)

    trajectory = json.loads(Path(trajectory_path).read_text(encoding="utf-8"))
    beta_a = arb(trajectory["segment"]["beta_a"])
    beta_b = arb(trajectory["segment"]["beta_b"])

    panel_states = PanelStates(trajectory, panel_index, ORDER)
    base = KoenigsBase(panel_states.beta, ORDER)

    # Hub modes: exact DFT of stored values, then polynomial interpolation in
    # the base.  No transcendental step, hence no remainder, hence the same
    # code path as the constant producer.
    beta_series = acb_series(
        [acb(panel_states.beta_mid), acb(panel_states.beta_half)], ORDER
    )
    base_series = 2 * (beta_series - beta_a) / (beta_b - beta_a) - 1
    hub_mode_models = {
        k: core.base_model(values, base_series)
        for k, values in core.build_hub_modes(trajectory).items()
    }

    rows = []
    unproved = []
    started = time.time()
    for x_index in range(node_count):
        x = arb(x_index) / RH_COUNT
        y, y_prime, y_second = panel_states.at(x)

        h = acb_series([], ORDER)
        h_prime = acb_series([], ORDER)
        h_second = acb_series([], ORDER)
        for k, mode in hub_mode_models.items():
            angle = 2 * arb.pi() * k * x
            exponential = acb(angle.cos(), angle.sin())
            h += mode * (exponential - 1)
            frequency = acb(0, 2 * arb.pi() * k)
            h_prime += mode * frequency * exponential
            h_second += mode * (-((2 * arb.pi() * k) ** 2)) * exponential
        h_beta = core.differentiate_beta_series(h, panel_states.beta_half)

        try:
            if x_index == 0:
                # T(0) = 1 is the removable endpoint of the Koenigs
                # coordinate.  Gamma(0) = 0 follows from the functional
                # equation Gamma(x) = c T(x) [Gamma(x-1) - T(x-1)] and
                # sexp_b(-1) = 0 (docs/UPPER_ENDPOINT_MODEL.md section 3);
                # Gamma'(0) = 0 is the stopped model's convention, kept here
                # so the derived and constant payloads differ only where a
                # remainder is charged.
                gamma = TM.zero(ORDER)
                gamma_prime = TM.zero(ORDER)
                koenigs_trace = {"endpoint_stopped_model": True, "depth": 0}
            else:
                gamma, gamma_prime, koenigs_trace = koenigs_pair(
                    base, y, state_radius, koenigs_stop, KOENIGS_MAX_DEPTH
                )
            tail, tail_y, tail_y_prime, tail_trace = certified_tail(
                base.c, y, y_prime, tail_stop
            )
        except (TaylorModelError, KoenigsError, TailError) as error:
            reason = f"{type(error).__name__}: {error}"
            unproved.append(
                {"panel": panel_index, "x_index": x_index, "reason": reason}
            )
            rows.append({"x_index": x_index, "unproved": reason})
            continue

        rows.append(
            {
                "x_index": x_index,
                "T": encode_model(y),
                "T_prime": encode_model(y_prime),
                "T_second": encode_model(y_second),
                "h": core.encode_series(h),
                "h_prime": core.encode_series(h_prime),
                "h_second": core.encode_series(h_second),
                "h_beta": core.encode_series(h_beta),
                "Gamma": encode_model(gamma),
                "Gamma_prime": encode_model(gamma_prime),
                "tail": encode_model(tail),
                "tail_partial_T": encode_model(tail_y),
                "tail_partial_T_prime": encode_model(tail_y_prime),
                "koenigs_trace": encode_trace(koenigs_trace),
                "tail_trace": encode_trace(tail_trace),
            }
        )
        if x_index and x_index % 24 == 0:
            print(
                f"[8r1-derived] panel {panel_index + 1}: x {x_index}/{node_count}"
                f"  {len(unproved)} unproved  {time.time() - started:.0f}s",
                flush=True,
            )

    return {
        "index": panel_index,
        "beta_left": core.outward_float_bounds(panel_states.beta_left),
        "beta_right": core.outward_float_bounds(panel_states.beta_right),
        "beta_midpoint": core.outward_float_bounds(panel_states.beta_mid),
        "beta_halfwidth": core.outward_float_bounds(panel_states.beta_half),
        "koenigs_base": {
            "psi_radius": core.outward_float_bounds(base.psi_radius),
            "psi_constant": core.outward_float_bounds(base.psi_constant),
            "phi_beta_constant": core.outward_float_bounds(
                base.phi_beta_constant
            ),
            "inverse_lambda_upper": core.outward_float_bounds(
                base.inverse_lambda_upper
            ),
            "fixed_point_remainder": core.outward_float_bounds(base.L.r),
        },
        "primitive_rows": rows,
        "unproved": unproved,
        "elapsed_seconds": round(time.time() - started, 1),
    }


def hull(records: list[list[str]]) -> list[str] | None:
    return core.interval_hull(records) if records else None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=PANEL_COUNT)
    parser.add_argument(
        "--nodes",
        type=int,
        default=RH_COUNT,
        help="nodes per panel (default all 96); smaller values for validation",
    )
    args = parser.parse_args()

    configure(ORDER, PRECISION_BITS)
    trajectory = json.loads(args.trajectory.read_text(encoding="utf-8"))
    if trajectory.get("schema_version") != "paper-vi-trajectory-8r1-v1":
        raise RuntimeError("unexpected trajectory schema")
    if len(trajectory.get("nodes", [])) != 7:
        raise RuntimeError("8r1 requires seven Lobatto states")
    if not 1 <= args.nodes <= RH_COUNT:
        raise RuntimeError(f"nodes must lie in 1..{RH_COUNT}")

    started = time.time()
    tasks = [
        (str(args.trajectory), panel_index, args.nodes)
        for panel_index in range(PANEL_COUNT)
    ]
    print(
        f"[8r1-derived] {PANEL_COUNT} panels x {args.nodes} nodes, "
        f"{args.jobs} processes",
        flush=True,
    )
    if args.jobs > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            panels = list(pool.map(produce_panel, tasks))
    else:
        panels = [produce_panel(task) for task in tasks]
    panels.sort(key=lambda panel: panel["index"])

    # Global gates are re-derived from the rows here, not taken from any
    # per-panel summary.  Endpoint rows run no orbit and are skipped; the
    # multiplier is a base property and comes from every panel.
    phase_margins: list[list[str]] = []
    terminal: list[list[str]] = []
    gamma_remainders: list[list[str]] = []
    gamma_prime_remainders: list[list[str]] = []
    gamma_model_remainders: list[list[str]] = []
    gamma_prime_model_remainders: list[list[str]] = []
    tail_remainders: list[list[str]] = []
    inverse_lambda: list[list[str]] = []
    unproved: list[dict] = []
    endpoint_rows = 0
    for panel in panels:
        inverse_lambda.append(panel["koenigs_base"]["inverse_lambda_upper"])
        unproved.extend(panel.pop("unproved"))
        for row in panel["primitive_rows"]:
            if "unproved" in row:
                continue
            trace = row["koenigs_trace"]
            if trace.get("endpoint_stopped_model"):
                endpoint_rows += 1
            else:
                phase_margins.append(trace["phase_margin_lower"])
                terminal.append(trace["terminal_distance"])
                gamma_remainders.append(trace["gamma_remainder"])
                gamma_prime_remainders.append(trace["gamma_prime_remainder"])
                gamma_model_remainders.append(trace["gamma_model_remainder"])
                gamma_prime_model_remainders.append(
                    trace["gamma_prime_model_remainder"]
                )
            tail_remainders.append(row["tail_trace"]["remainder"])

    source = Path(__file__).resolve()
    payload = {
        "schema_version": SCHEMA,
        "segment_id": "8r1",
        "status": "complete" if not unproved else "incomplete",
        "remainder_provenance": "derived (closure layer)",
        "interval_backend": {
            "name": "python-flint/Arb",
            "version": __import__("flint").__version__,
            "precision_bits": PRECISION_BITS,
            "series_cap": ORDER,
            "directed_rounding": True,
            "serialized_endpoints": "outward-widened IEEE-754 binary64",
        },
        "trajectory": {
            "path": args.trajectory.name,
            "file_sha256": core.sha256(args.trajectory),
            "payload_sha256": trajectory["payload_sha256"],
        },
        "conventions_sha256": trajectory["conventions_sha256"],
        "dimensions": {
            "base_lobatto_nodes": 7,
            "base_panels": PANEL_COUNT,
            "beta_series_order": ORDER,
            "state_chebyshev_count": STATE_COUNT,
            "rh_uniform_count": RH_COUNT,
            "rh_nodes_per_panel": args.nodes,
            "hub_mode_cutoff": core.HUB_MODE_CUTOFF,
        },
        "proof_constants": {
            "koenigs_stop": core.outward_float_bounds(arb(KOENIGS_STOP_TEXT)),
            "koenigs_max_depth": KOENIGS_MAX_DEPTH,
            "state_disc_radius": core.outward_float_bounds(
                arb(STATE_DISC_RADIUS_TEXT)
            ),
            "tail_stop": core.outward_float_bounds(arb(TAIL_STOP_TEXT)),
            "endpoint_model": (
                "Gamma(T(0))=0 from Gamma(x)=c*T(x)*[Gamma(x-1)-T(x-1)] and "
                "sexp_b(-1)=0; Gamma'(0)=0 is the stopped model's convention"
            ),
        },
        "derived_remainders": {
            "beta_taylor": (
                "certified Cauchy tail per elementary step; no set constant"
            ),
            "koenigs": (
                "telescoped psi/Phi_beta bounds on a certified disc "
                "(koenigs.koenigs_pair); Gamma' by a Cauchy estimate on the "
                "state disc; the factors 128 and 1024 are not used"
            ),
            "tail": (
                "proved ratio bound 1/(c y_{j+1}) with a checked "
                "forward-invariance premise (tails.certified_tail)"
            ),
            "state_and_hub": (
                "polynomial in the base and exact DFT of stored values; no "
                "transcendental step, no remainder in either producer"
            ),
        },
        "global_koenigs_gates": {
            "phase_margin_hull": hull(phase_margins),
            "inverse_lambda_hull": hull(inverse_lambda),
            "terminal_distance_hull": hull(terminal),
            "gamma_remainder_hull": hull(gamma_remainders),
            "gamma_prime_remainder_hull": hull(gamma_prime_remainders),
            "gamma_model_remainder_hull": hull(gamma_model_remainders),
            "gamma_prime_model_remainder_hull": hull(
                gamma_prime_model_remainders
            ),
            "tail_remainder_hull": hull(tail_remainders),
            "endpoint_rows": endpoint_rows,
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
                source.with_name("tails.py"),
                source.with_name("tmodel.py"),
                source.with_name("states.py"),
                source.with_name("certify_8r1.py"),
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

    covered = PANEL_COUNT * args.nodes
    print(f"[8r1-derived] wrote {args.output}", flush=True)
    print(f"[8r1-derived] payload sha256 {payload['payload_sha256']}")
    print(
        f"[8r1-derived] status {payload['status']}: {len(unproved)} of "
        f"{covered} rows unproved in {time.time() - started:.0f}s"
    )
    gates = payload["global_koenigs_gates"]
    for key in (
        "phase_margin_hull",
        "gamma_remainder_hull",
        "gamma_prime_remainder_hull",
        "tail_remainder_hull",
    ):
        print(f"[8r1-derived] {key:32s} {gates[key]}")
    return 0 if not unproved else 2


if __name__ == "__main__":
    raise SystemExit(main())
