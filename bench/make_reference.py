"""Generate frozen reference values with the ORIGINAL fatou.gp at dps+margin.

Run once after the 64-bit switch; afterwards research/reference/ is frozen.
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend import FatouGP, find_default_fatou_gp, find_default_gp_exe

from workload import all_cases, case_id

REF_MARGIN = 20
REFERENCE_PATH = Path(__file__).resolve().parents[1] / "research" / "reference" / "values.json"


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate frozen reference values")
    parser.add_argument("--deep", action="store_true", help="include dps 500/1000 tiers")
    parser.add_argument("--out", default=str(REFERENCE_PATH))
    args = parser.parse_args()

    fatou = find_default_fatou_gp()
    gp_exe = find_default_gp_exe()
    values: dict[str, dict[str, str]] = {}
    # group by (base, dps) so each sexpinit happens once per group
    groups: dict[tuple[str, int], list] = defaultdict(list)
    for c in all_cases(deep=args.deep):
        groups[(c.base, c.dps)].append(c)

    for (base, dps), cases in sorted(groups.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        ref_dps = dps + REF_MARGIN
        mp.mp.dps = ref_dps + 30
        gp = FatouGP(gp_exe=gp_exe, fatou_gp=fatou, dps=ref_dps,
                     looplim=max(35, ref_dps - 20))
        expressions = [f"{c.kind}({c.arg})" for c in cases]
        print(f"[make_reference] base={base} dps={dps} (ref_dps={ref_dps}) "
              f"{len(expressions)} exprs ...", flush=True)
        results = gp.eval_batch(base, expressions)
        for c, val in zip(cases, results):
            values[case_id(c)] = {
                "real": mp.nstr(mp.re(val), ref_dps, min_fixed=0, max_fixed=0),
                "imag": mp.nstr(mp.im(val), ref_dps, min_fixed=0, max_fixed=0),
            }

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    existing: dict[str, dict[str, str]] = {}
    if out.exists():
        existing = json.loads(out.read_text(encoding="utf-8")).get("values", {})
    existing.update(values)
    payload = {
        "meta": {
            "generated": datetime.date.today().isoformat(),
            "gp_exe": str(gp_exe),
            "fatou_sha256": hashlib.sha256(Path(fatou).read_bytes()).hexdigest(),
            "ref_dps_margin": REF_MARGIN,
        },
        "values": existing,
    }
    out.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"[make_reference] wrote {len(values)} values to {out}")


if __name__ == "__main__":
    main()
