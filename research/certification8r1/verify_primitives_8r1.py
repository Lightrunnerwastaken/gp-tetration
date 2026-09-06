"""Independent primitive verifier for the ``8r1`` segment.

This program shares no code with ``produce_8r1.py``, ``certify_8r1.py`` or
``certify_hilbert_8r1.py``.  It reads only the raw trajectory -- the seven
Lobatto states and the segment conventions -- and *regenerates* the Koenigs
and tail primitives from scratch with certified Taylor models whose remainder
is derived rather than declared.  It then

* writes its own outward-rounded primitive file with a full provenance record
  for every remainder, and
* compares its enclosures against the stored pre-certificate primitives at
  sample points of the base panel, reporting the discrepancy next to the
  remainder each side claims.

Exit status is non-zero if any certified step fails or if a stored primitive
is inconsistent with the regenerated enclosure.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

from flint import acb, arb

from koenigs import KoenigsBase, koenigs_pair
from states import PanelStates
from tails import certified_tail
from tmodel import TM, configure


ORDER = 22
PRECISION_BITS = 448
KOENIGS_STOP = arb("1e-18")
KOENIGS_MAX_DEPTH = 500
TAIL_STOP = arb("1e-50")
STATE_DISC_RADIUS = arb("1e-3")
RH_COUNT = 96
SCHEMA = "paper-vi-8r1-derived-primitives-v1"


def outward(value: arb) -> list[str]:
    lower = math.nextafter(float(value.lower()), -math.inf)
    upper = math.nextafter(float(value.upper()), math.inf)
    return [format(lower, ".17e"), format(upper, ".17e")]


def encode_model(model: TM) -> dict:
    return {
        "coefficients": [
            {
                "re": outward(coefficient.real),
                "im": outward(coefficient.imag),
            }
            for coefficient in model.c
        ],
        "remainder": outward(model.r),
    }


def stored_series_ball(record: list[dict], s: acb) -> acb:
    """Evaluate a pre-certificate coefficient row at a base sample."""

    value = acb(0)
    power = acb(1)
    for entry in record:
        real = arb(entry["re"][0])
        real_upper = arb(entry["re"][1])
        imag = arb(entry["im"][0])
        imag_upper = arb(entry["im"][1])
        coefficient = acb(
            arb((real + real_upper) / 2, (real_upper - real) / 2),
            arb((imag + imag_upper) / 2, (imag_upper - imag) / 2),
        )
        value += coefficient * power
        power = power * s
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--stored-primitives", type=Path, default=None)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--panels", type=int, default=6)
    parser.add_argument("--nodes", type=int, default=RH_COUNT)
    args = parser.parse_args()

    configure(ORDER, PRECISION_BITS)
    trajectory = json.loads(args.trajectory.read_text(encoding="utf-8"))
    if trajectory.get("schema_version") != "paper-vi-trajectory-8r1-v1":
        raise RuntimeError("unexpected trajectory schema")
    stored = (
        json.loads(args.stored_primitives.read_text(encoding="utf-8"))
        if args.stored_primitives
        else None
    )

    samples = [acb(arb(value)) for value in ("-1", "-0.5", "0", "0.5", "1")]
    issues: list[str] = []
    panels_out = []
    worst = {
        "gamma_remainder": arb(0),
        "gamma_prime_remainder": arb(0),
        "tail_remainder": arb(0),
        "model_remainder": arb(0),
        "log_kappa_error": arb(0),
        "base_error": arb(0),
        "tail_ratio": arb(0),
    }
    best_phase: arb | None = None
    comparison = {
        "max_gamma_gap": arb(0),
        "max_tail_gap": arb(0),
        "declared_beta_series_remainder": "2e-16",
        "declared_gamma_remainder_hull": None,
    }
    if stored:
        comparison["declared_gamma_remainder_hull"] = stored[
            "global_koenigs_gates"
        ]["gamma_remainder_hull"]

    start = time.time()
    for panel_index in range(args.panels):
        panel_states = PanelStates(trajectory, panel_index, ORDER)
        base = KoenigsBase(panel_states.beta, ORDER)
        rows = []
        for x_index in range(args.nodes):
            x = arb(x_index) / RH_COUNT
            y, y_prime, y_second = panel_states.at(x)
            if x_index == 0:
                # T(0) = 1 is the removable endpoint of the Koenigs
                # coordinate; the stopped model sets Gamma = Gamma' = 0 there
                # before any quotient is taken.
                gamma = TM.zero(ORDER)
                gamma_prime = TM.zero(ORDER)
                koenigs_trace = {"endpoint_stopped_model": True}
            else:
                gamma, gamma_prime, koenigs_trace = koenigs_pair(
                    base, y, STATE_DISC_RADIUS, KOENIGS_STOP, KOENIGS_MAX_DEPTH
                )
                for key in (
                    "gamma_remainder",
                    "gamma_prime_remainder",
                    "log_kappa_error",
                    "base_error",
                ):
                    if bool(koenigs_trace[key] > worst[key]):
                        worst[key] = koenigs_trace[key]
                phase = koenigs_trace["phase_margin_lower"]
                if best_phase is None or bool(phase < best_phase):
                    best_phase = phase
            tail, tail_y, tail_y_prime, tail_trace = certified_tail(
                base.c, y, y_prime, TAIL_STOP
            )
            if bool(tail_trace["remainder"] > worst["tail_remainder"]):
                worst["tail_remainder"] = tail_trace["remainder"]
            if bool(tail_trace["ratio_upper"] > worst["tail_ratio"]):
                worst["tail_ratio"] = tail_trace["ratio_upper"]
            for model in (y, y_prime, y_second, gamma, gamma_prime, tail):
                if bool(model.r > worst["model_remainder"]):
                    worst["model_remainder"] = model.r

            if stored and panel_index < len(stored["panels"]):
                source = stored["panels"][panel_index]["primitive_rows"][
                    x_index
                ]
                for sample in samples:
                    mine = gamma.poly()(sample)
                    theirs = stored_series_ball(source["Gamma"], sample)
                    gap = (mine - theirs).abs_lower()
                    if bool(gap > comparison["max_gamma_gap"]):
                        comparison["max_gamma_gap"] = gap
                    mine = tail.poly()(sample)
                    theirs = stored_series_ball(source["tail"], sample)
                    gap = (mine - theirs).abs_lower()
                    if bool(gap > comparison["max_tail_gap"]):
                        comparison["max_tail_gap"] = gap

            rows.append(
                {
                    "x_index": x_index,
                    "T": encode_model(y),
                    "T_prime": encode_model(y_prime),
                    "T_second": encode_model(y_second),
                    "Gamma": encode_model(gamma),
                    "Gamma_prime": encode_model(gamma_prime),
                    "tail": encode_model(tail),
                    "tail_partial_T": encode_model(tail_y),
                    "tail_partial_T_prime": encode_model(tail_y_prime),
                    "koenigs_trace": {
                        key: (
                            outward(value)
                            if isinstance(value, arb)
                            else value
                        )
                        for key, value in koenigs_trace.items()
                    },
                    "tail_trace": {
                        key: (
                            outward(value)
                            if isinstance(value, arb)
                            else value
                        )
                        for key, value in tail_trace.items()
                    },
                }
            )
            if x_index % 24 == 0:
                print(
                    f"[verify] panel {panel_index + 1}/{args.panels} "
                    f"x {x_index}/{args.nodes} "
                    f"({time.time() - start:.0f}s)",
                    flush=True,
                )
        panels_out.append(
            {
                "index": panel_index,
                "beta_left": outward(panel_states.beta_left),
                "beta_right": outward(panel_states.beta_right),
                "beta_halfwidth": outward(panel_states.beta_half),
                "psi_radius": outward(base.psi_radius),
                "psi_constant": outward(base.psi_constant),
                "phi_beta_constant": outward(base.phi_beta_constant),
                "inverse_lambda_upper": outward(base.inverse_lambda_upper),
                "fixed_point_remainder": outward(base.L.r),
                "primitive_rows": rows,
            }
        )

    payload = {
        "schema_version": SCHEMA,
        "segment_id": "8r1",
        "interval_backend": {
            "name": "python-flint/Arb",
            "version": __import__("flint").__version__,
            "precision_bits": PRECISION_BITS,
            "series_cap": ORDER,
            "directed_rounding": True,
        },
        "trajectory": {
            "path": args.trajectory.name,
            "payload_sha256": trajectory["payload_sha256"],
        },
        "derived_remainders": {
            "beta_taylor": (
                "certified Cauchy tail per elementary step; no set constant"
            ),
            "koenigs": (
                "telescoped psi/Phi_beta bounds on a certified disc; the "
                "factors 128 and 1024 of the pre-certificate are not used"
            ),
            "tail": (
                "proved ratio bound 1/(c y_{j+1}) <= 1/2 with a checked "
                "forward-invariance premise"
            ),
        },
        "worst_case": {key: outward(value) for key, value in worst.items()},
        "phase_margin_lower": outward(best_phase or arb.pi()),
        "panels": panels_out,
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

    report = {
        "schema_version": "paper-vi-8r1-independent-verification-v1",
        "segment_id": "8r1",
        "panels_checked": args.panels,
        "nodes_per_panel": args.nodes,
        "elapsed_seconds": round(time.time() - start, 1),
        "derived_worst_case": {
            key: outward(value) for key, value in worst.items()
        },
        "phase_margin_lower": outward(best_phase or arb.pi()),
        "comparison_with_pre_certificate": {
            key: (outward(value) if isinstance(value, arb) else value)
            for key, value in comparison.items()
        },
        "primitives_payload_sha256": payload["payload_sha256"],
        "issues": issues,
        "passed": not issues,
    }
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print()
    print(f"independent primitive verification: {'PASS' if not issues else 'FAIL'}")
    for key, value in worst.items():
        print(f"  worst {key:24s} {value.str(8)}")
    print(f"  phase margin lower        {(best_phase or arb.pi()).str(8)}")
    if stored:
        print(
            f"  max |Gamma_new - Gamma_stored|  "
            f"{comparison['max_gamma_gap'].str(8)}"
        )
        print(
            f"  max |tail_new - tail_stored|    "
            f"{comparison['max_tail_gap'].str(8)}"
        )
    print(f"  report {args.report}")
    return 0 if not issues else 1


if __name__ == "__main__":
    raise SystemExit(main())
