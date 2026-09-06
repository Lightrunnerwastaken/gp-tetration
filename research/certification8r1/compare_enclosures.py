"""Containment check between the derived and the stored ``8r1`` primitives.

The independent verifier writes coefficient rows *together with the derived
remainder* of each model.  The pre-certificate writes coefficient rows only:
its single set constant ``BETA_SERIES_REMAINDER`` was folded into the
constant coefficient at production time.  A meaningful comparison is
therefore not "are the coefficients equal" but "do the two enclosures of the
same function overlap".

For each model and each base sample ``s`` this program forms

    gap(s) = min over the two coefficient balls of | P_derived(s) - P_stored(s) |

and reports it against the derived remainder.  ``gap <= remainder`` means the
stored enclosure meets the certified one, so the pre-certificate's primitive
is consistent with the proved bound.  A gap larger than the remainder would
mean the two disagree by more than the certified tolerance allows -- that is
the failure the check exists to catch.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from flint import acb, arb, ctx


KEYS = ("T", "T_prime", "T_second", "Gamma", "Gamma_prime", "tail")


def ball_from_pair(pair: list[str]) -> arb:
    lower = arb(pair[0])
    upper = arb(pair[1])
    return arb((lower + upper) / 2, (upper - lower) / 2)


def evaluate_derived(model: dict, s: acb) -> tuple[acb, arb]:
    value = acb(0)
    power = acb(1)
    for entry in model["coefficients"]:
        coefficient = acb(
            ball_from_pair(entry["re"]), ball_from_pair(entry["im"])
        )
        value += coefficient * power
        power = power * s
    return value, ball_from_pair(model["remainder"]).upper()


def evaluate_stored(record: list[dict], s: acb) -> acb:
    value = acb(0)
    power = acb(1)
    for entry in record:
        coefficient = acb(
            ball_from_pair(entry["re"]), ball_from_pair(entry["im"])
        )
        value += coefficient * power
        power = power * s
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--derived", type=Path, required=True)
    parser.add_argument("--stored", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    ctx.prec = 320
    derived = json.loads(args.derived.read_text(encoding="utf-8"))
    stored = json.loads(args.stored.read_text(encoding="utf-8"))
    samples = [acb(arb(value)) for value in ("-1", "-0.5", "0", "0.5", "1")]

    summary = {}
    failures = []
    for key in KEYS:
        worst_gap = arb(0)
        worst_remainder = arb(0)
        worst_location = None
        for panel_derived, panel_stored in zip(
            derived["panels"], stored["panels"]
        ):
            for row_derived, row_stored in zip(
                panel_derived["primitive_rows"], panel_stored["primitive_rows"]
            ):
                if key not in row_stored:
                    continue
                for index, s in enumerate(samples):
                    mine, remainder = evaluate_derived(row_derived[key], s)
                    theirs = evaluate_stored(row_stored[key], s)
                    gap = (mine - theirs).abs_lower()
                    if bool(gap > worst_gap):
                        worst_gap = gap
                        worst_remainder = remainder
                        worst_location = {
                            "panel": panel_derived["index"],
                            "x_index": row_derived["x_index"],
                            "sample": index,
                        }
                    if not bool(gap <= remainder):
                        failures.append(
                            {
                                "quantity": key,
                                "panel": panel_derived["index"],
                                "x_index": row_derived["x_index"],
                                "sample": index,
                                "gap": gap.str(12),
                                "derived_remainder": remainder.str(12),
                            }
                        )
        summary[key] = {
            "worst_gap": worst_gap.str(12),
            "derived_remainder_there": worst_remainder.str(12),
            "contained": bool(worst_gap <= worst_remainder),
            "location": worst_location,
        }
        status = "ok  " if summary[key]["contained"] else "FAIL"
        print(
            f"  [{status}] {key:12s} worst gap {worst_gap.str(8):>18s}  "
            f"derived remainder {worst_remainder.str(8)}"
        )

    report = {
        "schema_version": "paper-vi-8r1-enclosure-containment-v1",
        "segment_id": "8r1",
        "derived_payload_sha256": derived.get("payload_sha256"),
        "stored_payload_sha256": stored.get("payload_sha256"),
        "samples": ["-1", "-0.5", "0", "0.5", "1"],
        "per_quantity": summary,
        "failures": failures[:64],
        "failure_count": len(failures),
        "passed": not failures,
    }
    args.report.write_text(
        json.dumps(report, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print()
    print(
        f"enclosure containment: {'PASS' if not failures else 'FAIL'} "
        f"({len(failures)} violations)"
    )
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
