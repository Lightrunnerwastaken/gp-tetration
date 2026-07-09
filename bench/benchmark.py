"""Benchmark: verified correct digits per second on the frozen workload.

cold  = fresh FatouGP per group, timing includes process start + sexpinit
warm  = same FatouGP object, second evaluation of the group (persistent
        workers make this meaningful; on the one-shot path warm equals cold
        and is reported anyway for comparability)
"""
from __future__ import annotations

import argparse
import datetime
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend import FatouGP, find_default_fatou_gp, find_default_gp_exe

from metrics import correct_digits, load_reference
from workload import all_cases, case_id

REFERENCE_PATH = Path(__file__).resolve().parents[1] / "research" / "reference" / "values.json"
RESULTS_DIR = Path(__file__).resolve().parent / "results"


def run_group(gp: FatouGP, base: str, cases: list, reference: dict) -> dict:
    expressions = [f"{c.kind}({c.arg})" for c in cases]
    start = time.perf_counter()
    results = gp.eval_batch(base, expressions)
    elapsed = time.perf_counter() - start
    dps = cases[0].dps
    digits = [correct_digits(v, reference[case_id(c)], cap=dps - 10)
              for c, v in zip(cases, results)]
    min_digits = min(digits)
    return {
        "seconds": elapsed,
        "values": len(cases),
        "min_correct_digits": min_digits,
        "digits_per_second": sum(digits) / elapsed,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="fatou_backend benchmark")
    parser.add_argument("--mode", choices=["all", "throughput", "precision"], default="all")
    parser.add_argument("--repeat", type=int, default=1)
    parser.add_argument("--deep", action="store_true")
    parser.add_argument("--fatou", default=None, help="fatou.gp variant to benchmark")
    parser.add_argument("--label", default="run")
    parser.add_argument("--json", default=None)
    args = parser.parse_args()

    fatou = Path(args.fatou) if args.fatou else find_default_fatou_gp()
    gp_exe = find_default_gp_exe()
    if not REFERENCE_PATH.exists():
        raise SystemExit(f"reference file missing: {REFERENCE_PATH}")

    tiers = {"all": {"throughput", "precision"}, "throughput": {"throughput"},
             "precision": {"precision"}}[args.mode]
    groups: dict[tuple[str, int], list] = defaultdict(list)
    for c in all_cases(deep=args.deep):
        if c.tier in tiers:
            groups[(c.base, c.dps)].append(c)

    runs = []
    for (base, dps), cases in sorted(groups.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        mp.mp.dps = dps + 40
        reference = load_reference(REFERENCE_PATH)  # reparse at sufficient dps
        for rep in range(args.repeat):
            gp = FatouGP(gp_exe=gp_exe, fatou_gp=fatou, dps=dps,
                         looplim=max(35, dps - 20))
            tier = cases[0].tier
            cold = run_group(gp, base, cases, reference)
            cold.update({"group": f"{tier}|{base}|{dps}", "mode": "cold", "rep": rep})
            runs.append(cold)
            warm = run_group(gp, base, cases, reference)
            warm.update({"group": f"{tier}|{base}|{dps}", "mode": "warm", "rep": rep})
            runs.append(warm)
            close = getattr(gp, "close", None)
            if close is not None:
                close()
            print(f"{tier}|{base}|{dps} rep {rep}: "
                  f"cold {cold['seconds']:.2f}s ({cold['digits_per_second']:.0f} d/s), "
                  f"warm {warm['seconds']:.2f}s ({warm['digits_per_second']:.0f} d/s), "
                  f"min digits {min(cold['min_correct_digits'], warm['min_correct_digits']):.1f}",
                  flush=True)

    payload = {
        "meta": {"label": args.label, "date": datetime.datetime.now().isoformat(timespec="seconds"),
                 "gp_exe": str(gp_exe), "fatou": str(fatou), "repeat": args.repeat},
        "runs": runs,
    }
    out = Path(args.json) if args.json else RESULTS_DIR / f"{datetime.date.today().isoformat()}-{args.label}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"[benchmark] wrote {out}")


if __name__ == "__main__":
    main()
