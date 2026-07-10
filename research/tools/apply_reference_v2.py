"""Apply the approved v2 reference correction (HUMAN-APPROVED runs only).

Run with --confirm "Option A" after the user approved
research/DECISION_reference_v2.md. Steps:
  1. Replace base-e entries in research/reference/values.json with the
     verified values from values_v2_proposed.json (provenance recorded).
  2. Recalibrate gate REQUIRED_DIGITS for base-e groups: measure the current
     fork against the corrected references, set required = worst - 5.
  3. Print the gate.py edits to make (the script does NOT edit gate.py;
     that final step stays manual/reviewed, then the file is re-frozen).
"""
from __future__ import annotations

import argparse
import datetime
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

VALUES = REPO / "research" / "reference" / "values.json"
PROPOSED = REPO / "research" / "reference" / "values_v2_proposed.json"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--confirm", required=True,
                        help='must be exactly "Option A" (user approval token)')
    args = parser.parse_args()
    if args.confirm != "Option A":
        raise SystemExit("refusing: pass --confirm \"Option A\" after user approval")

    proposed = json.loads(PROPOSED.read_text(encoding="utf-8"))
    payload = json.loads(VALUES.read_text(encoding="utf-8"))
    replaced = []
    for key, entry in proposed["values"].items():
        if key in payload["values"]:
            payload["values"][key] = entry
            replaced.append(key)
    payload["meta"]["v2_correction"] = {
        "date": datetime.date.today().isoformat(),
        "replaced": sorted(replaced),
        "verification_digits": proposed.get("verification_digits", {}),
        "decision": "research/DECISION_reference_v2.md Option A",
    }
    VALUES.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"[apply_v2] replaced {len(replaced)} base-e entries in values.json")

    # step 2: measure current fork against corrected references
    sys.path.insert(0, str(REPO / "research"))
    from gate import measure_agreement, DEFAULT_FORK
    agreement, errors = measure_agreement(DEFAULT_FORK)
    for line in errors:
        print(line)
    print("[apply_v2] gate.py REQUIRED_DIGITS updates for base-e groups:")
    for group, worst in sorted(agreement.items()):
        if group.startswith("e|"):
            print(f'    "{group}": {max(1.0, worst - 5):.1f},   # was-measured {worst:.1f}')
    print("[apply_v2] edit gate.py accordingly, run the slow gate tests, commit.")


if __name__ == "__main__":
    main()
