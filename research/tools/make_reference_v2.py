"""Generate PROPOSED corrected reference values for base-e cases.

Root cause (journal 2026-07-10): the frozen references were produced with
nlim=30; base e converges at ~2.12 true digits/iteration, so every base-e
reference is only ~64 digits accurate regardless of its nominal dps.

This script recomputes base-e cases with a converged engine (nlim scaled to
the digits goal, looplim=0 = engine-internal full target) and VERIFIES each
value by decorrelation: an independent run with ~1.5x nlim at dps+20 must
agree beyond the digits needed. Output goes to a NEW file
(values_v2_proposed.json) — research/reference/values.json stays frozen until
the human decision.
"""
from __future__ import annotations

import datetime
import json
import sys
import time
from pathlib import Path

import mpmath as mp

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

from bench.metrics import correct_digits
from bench.workload import all_cases, case_id
from fatou_backend import FatouGP, find_default_fatou_gp, find_default_gp_exe

OUT_PATH = REPO / "research" / "reference" / "values_v2_proposed.json"
RATE_E = 2.1  # measured true digits per iteration for base e


def converged_run(dps: int, nlim: int, expressions: list[str],
                  base: str = "e") -> list[mp.mpc]:
    gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp=find_default_fatou_gp(),
                 dps=dps, nlim=nlim, nskip=4, looplim=0, init_timeout=4 * 3600.0)
    try:
        return gp.eval_batch(base, expressions)
    finally:
        gp.close()


def main() -> None:
    # (base, case-dps filter, ref_dps, nlim main, verify dps, verify nlim)
    tiers = [
        ("e", 80, 100, 60, 120, 90),
        # retry 2: the loop self-terminates at contour-re ~ precis-throwp;
        # TRUE accuracy lags contour-re by a scale-dependent gap (~27 at
        # dps 220). nlim>=~130 is non-binding; raise WORKING dps instead.
        ("e", 200, 250, 220, 270, 260),
    ]
    if "--dps500" in sys.argv:
        # corrected model: loop exits at contour ~ precis-throwp; truth lags
        # by a scale gap -> main working dps 560 for >=510 true digits
        tiers = [("e", 500, 560, 340, 580, 360)]
    if "--base2-1000" in sys.argv:
        # base 2 converges ~32 digits/iter; nlim=30 capped it at ~956 true
        # digits, the frozen 2|1000 reference claims 1020 -> regenerate with
        # non-binding nlim and raised working dps
        tiers = [("2", 1000, 1060, 60, 1080, 70)]
    values: dict[str, dict[str, str]] = {}
    verification: dict[str, float] = {}
    existing = {}
    if OUT_PATH.exists():
        payload = json.loads(OUT_PATH.read_text(encoding="utf-8"))
        existing = payload.get("values", {})
        verification = payload.get("verification_digits", {})

    for base, case_dps, ref_dps, nlim_main, ver_dps, ver_nlim in tiers:
        cases = [c for c in all_cases(deep=True) if c.base == base and c.dps == case_dps]
        expressions = [f"{c.kind}({c.arg})" for c in cases]
        needed = case_dps + 10  # ref must be true to case_dps + margin
        mp.mp.dps = ver_dps + 60
        t = time.perf_counter()
        print(f"[v2] {base}|{case_dps}: main run dps={ref_dps} nlim={nlim_main} "
              f"({len(cases)} exprs) ...", flush=True)
        main_vals = converged_run(ref_dps, nlim_main, expressions, base=base)
        print(f"[v2]   main done ({time.perf_counter()-t:.0f}s); verify run "
              f"dps={ver_dps} nlim={ver_nlim} ...", flush=True)
        t = time.perf_counter()
        ver_vals = converged_run(ver_dps, ver_nlim, expressions, base=base)
        print(f"[v2]   verify done ({time.perf_counter()-t:.0f}s)", flush=True)
        worst = min(correct_digits(m, v, cap=ref_dps)
                    for m, v in zip(main_vals, ver_vals))
        status = "OK" if worst >= needed else "INSUFFICIENT"
        print(f"[v2] {base}|{case_dps}: verified to {worst:.1f} digits "
              f"(needed {needed}) -> {status}", flush=True)
        if worst < needed:
            continue
        for c, val in zip(cases, main_vals):
            values[case_id(c)] = {
                "real": mp.nstr(mp.re(val), ref_dps, min_fixed=0, max_fixed=0),
                "imag": mp.nstr(mp.im(val), ref_dps, min_fixed=0, max_fixed=0),
            }
            verification[case_id(c)] = round(worst, 1)

    existing.update(values)
    OUT_PATH.write_text(json.dumps({
        "meta": {
            "generated": datetime.date.today().isoformat(),
            "note": "PROPOSED corrected base-e references (converged nlim); "
                    "pending human approval to replace frozen values",
            "rate_digits_per_iter_base_e": RATE_E,
        },
        "values": existing,
        "verification_digits": verification,
    }, indent=1), encoding="utf-8")
    print(f"[v2] wrote {len(values)} new values (total {len(existing)}) to {OUT_PATH}")


if __name__ == "__main__":
    main()
