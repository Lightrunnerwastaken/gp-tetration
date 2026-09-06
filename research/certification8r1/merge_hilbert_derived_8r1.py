"""Merge sharded ``certify_hilbert_derived_8r1.py`` runs into one payload.

The derived producer takes ``--panel-start`` / ``--node-start``, so a full
chain can be run as six concurrent panel shards instead of one 2.5-hour
sequence.  This tool reassembles them and re-derives the global gates from the
rows, rather than trusting each shard's own summary.

Refusals are the point of the derived producer, so they are merged too: a
merged payload is ``complete`` only if every shard covered its full node range
and no row was refused.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from flint import arb

import certify_8r1 as core


SCHEMA = "paper-vi-8r1-hilbert-primitives-derived-v1"


def lower(bounds: list[str]) -> arb:
    return arb(bounds[0])


def upper(bounds: list[str]) -> arb:
    return arb(bounds[1])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shard", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expect-panels", type=int, default=6)
    parser.add_argument("--expect-nodes", type=int, default=864)
    args = parser.parse_args()

    shards = [
        json.loads(path.read_text(encoding="utf-8")) for path in args.shard
    ]
    for path, shard in zip(args.shard, shards):
        if shard.get("schema_version") != SCHEMA:
            raise RuntimeError(f"{path.name}: unexpected schema")
        recomputed = hashlib.sha256(
            json.dumps(
                {k: v for k, v in shard.items() if k != "payload_sha256"},
                sort_keys=True,
                separators=(",", ":"),
                ensure_ascii=False,
            ).encode("utf-8")
        ).hexdigest()
        if recomputed != shard.get("payload_sha256"):
            raise RuntimeError(f"{path.name}: payload hash mismatch")

    trajectory_hashes = {
        shard["trajectory"]["file_sha256"] for shard in shards
    }
    if len(trajectory_hashes) != 1:
        raise RuntimeError("shards were produced from different trajectories")

    panels: dict[int, dict] = {}
    for shard in shards:
        for panel in shard["panels"]:
            if panel["index"] in panels:
                raise RuntimeError(f"panel {panel['index']} appears twice")
            panels[panel["index"]] = panel

    missing = sorted(set(range(args.expect_panels)) - set(panels))
    if missing:
        raise RuntimeError(f"panels missing from the shard set: {missing}")

    # Re-derive the gates and the refusal list from the rows themselves.  A
    # shard summary that disagrees with its own rows would otherwise survive
    # the merge.
    unproved: list[dict] = []
    short_panels: list[int] = []
    phase_margin: arb | None = None
    inverse_lambda = arb(0)
    terminal = arb(0)
    endpoint_rows = 0
    endpoint_charge = arb(0)

    for index in sorted(panels):
        rows = panels[index]["trace_rows"]
        if len(rows) != args.expect_nodes:
            short_panels.append(index)
        for row in rows:
            if row.get("u") is None:
                unproved.append(
                    {
                        "panel": index,
                        "quadrature_index": row["quadrature_index"],
                        "reason": row.get("unproved", "unknown"),
                    }
                )
                continue
            trace = row["koenigs_trace"]
            # Endpoint rows run no orbit; they have no phase margin and no
            # multiplier to fold into the gates, only their own charge.
            if trace.get("endpoint_stopped_model"):
                endpoint_rows += 1
                value = upper(trace["endpoint_charge"])
                if bool(value > endpoint_charge):
                    endpoint_charge = value
                continue
            value = lower(trace["phase_margin_lower"])
            if phase_margin is None or bool(value < phase_margin):
                phase_margin = value
            value = upper(trace["inverse_lambda_upper"])
            if bool(value > inverse_lambda):
                inverse_lambda = value
            value = upper(trace["terminal_distance"])
            if bool(value > terminal):
                terminal = value

    reference = shards[0]
    payload = {
        "schema_version": SCHEMA,
        "segment_id": "8r1",
        "status": (
            "complete" if not unproved and not short_panels else "incomplete"
        ),
        "remainder_provenance": "derived (closure layer)",
        "merged_from_shards": len(shards),
        "interval_backend": reference["interval_backend"],
        "trajectory": reference["trajectory"],
        "dimensions": {
            "base_panels": len(panels),
            "beta_series_order": reference["dimensions"]["beta_series_order"],
            "quadrature_nodes": args.expect_nodes,
            "quadrature_order": reference["dimensions"]["quadrature_order"],
        },
        "koenigs_settings": reference["koenigs_settings"],
        "upper_endpoint_model": {
            **reference["upper_endpoint_model"],
            "rows": endpoint_rows,
            "worst_charge": core.outward_float_bounds(endpoint_charge),
        },
        "quadrature_proof": reference["quadrature_proof"],
        "global_gates": {
            "phase_margin_lower": core.outward_float_bounds(
                phase_margin if phase_margin is not None else arb(0)
            ),
            "inverse_lambda_upper": core.outward_float_bounds(inverse_lambda),
            "terminal_distance_upper": core.outward_float_bounds(terminal),
        },
        "short_panels": short_panels,
        "unproved_rows": unproved,
        "unproved_row_count": len(unproved),
        "panels": [panels[index] for index in sorted(panels)],
        "source_manifest": reference["source_manifest"],
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

    total = len(panels) * args.expect_nodes
    print(f"[merge] wrote {args.output}")
    print(f"[merge] payload sha256 {payload['payload_sha256']}")
    print(
        f"[merge] status {payload['status']}: "
        f"{len(unproved)} of {total} rows unproved"
    )
    print(
        f"[merge] upper endpoint model on {endpoint_rows} rows, "
        f"worst charge {float(endpoint_charge):.3e}"
    )
    if short_panels:
        print(f"[merge] panels with an incomplete node range: {short_panels}")
    return 0 if payload["status"] == "complete" else 2


if __name__ == "__main__":
    raise SystemExit(main())
