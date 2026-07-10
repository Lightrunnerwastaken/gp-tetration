"""v4 reference sanitation: regenerate ALL fast-base cases with converged
runs (looplim=0, nlim non-binding) and verify each group by the ERROR-VECTOR
method — the only proven verification (journal 2026-07-10): two runs with
different iteration budgets must agree beyond the needed digits; their delta
equals the worse run's true error (sexp error == contour error, factor 1).

Writes directly into research/reference/values.json (delegated authority),
recording provenance in meta.v4_correction.
"""
from __future__ import annotations

import datetime
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import mpmath as mp

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

from bench.metrics import correct_digits
from bench.workload import all_cases, case_id
from fatou_backend import FatouGP, find_default_fatou_gp, find_default_gp_exe

VALUES = REPO / "research" / "reference" / "values.json"
FAST_BASES = {"2", "10", "1+I", "0.8+0.4*I", "2+I"}


def converged_run(base: str, dps: int, nlim: int, expressions: list[str]) -> list[mp.mpc]:
    gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp=find_default_fatou_gp(),
                 dps=dps, nlim=nlim, nskip=4, looplim=0, init_timeout=4 * 3600.0)
    try:
        return gp.eval_batch(base, expressions)
    finally:
        gp.close()


def main() -> None:
    groups: dict[tuple[str, int], list] = defaultdict(list)
    for c in all_cases(deep=False):
        if c.base in FAST_BASES:
            groups[(c.base, c.dps)].append(c)

    payload = json.loads(VALUES.read_text(encoding="utf-8"))
    prov = payload["meta"].setdefault("v4_correction", {})

    for (base, case_dps), cases in sorted(groups.items(), key=lambda kv: (kv[0][1], kv[0][0])):
        ref_dps = case_dps + 20
        needed = case_dps + 10
        # both runs fully converge (loop exits at looplim ~ their dps), but
        # to DIFFERENT iteration depths: the delta equals the worse run's
        # true error (~ref_dps level) — proving >= needed digits without
        # sharing an nlim-cap systematic
        nlim_gen = 400
        mp.mp.dps = ref_dps + 80
        expressions = [f"{c.kind}({c.arg})" for c in cases]
        t = time.perf_counter()
        print(f"[v4] {base}|{case_dps}: main (dps {ref_dps}) ...", flush=True)
        main_vals = converged_run(base, ref_dps, nlim_gen, expressions)
        print(f"[v4]   main {time.perf_counter()-t:.0f}s; check (dps {ref_dps+13}) ...", flush=True)
        t = time.perf_counter()
        check_vals = converged_run(base, ref_dps + 13, nlim_gen, expressions)
        worst = min(correct_digits(m, v, cap=ref_dps)
                    for m, v in zip(main_vals, check_vals))
        status = "OK" if worst >= needed else "INSUFFICIENT"
        print(f"[v4] {base}|{case_dps}: fehler-vektor {worst:.1f} digits "
              f"(needed {needed}, check {time.perf_counter()-t:.0f}s) -> {status}", flush=True)
        if worst < needed:
            continue
        for c, val in zip(cases, main_vals):
            payload["values"][case_id(c)] = {
                "real": mp.nstr(mp.re(val), ref_dps, min_fixed=0, max_fixed=0),
                "imag": mp.nstr(mp.im(val), ref_dps, min_fixed=0, max_fixed=0),
            }
        prov[f"{base}|{case_dps}"] = {
            "date": datetime.date.today().isoformat(),
            "method": f"Original looplim=0 nlim={nlim_main} dps={ref_dps}; "
                      f"Fehler-Vektor vs nlim={nlim_check}: {worst:.1f} digits",
        }
        VALUES.write_text(json.dumps(payload, indent=1), encoding="utf-8")
        print(f"[v4]   {len(cases)} Werte geschrieben", flush=True)

    print("[v4] fertig")


if __name__ == "__main__":
    main()
