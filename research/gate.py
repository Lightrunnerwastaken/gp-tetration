"""Immutable correctness gate for fatou_fork.gp experiments.

FROZEN after task 9 of the 2026-07-09 plan. Loop experiments must never edit
this file or research/reference/. If the gate itself needs a change, that is a
human decision outside the loop, documented in the journal.

Thresholds were calibrated ONCE against the unmodified original engine
(python research/gate.py --calibrate) and then frozen: required digits per
group = original worst agreement - 5 safety digits; roundtrip max = original
worst residual * 100.
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import mpmath as mp

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

from bench.metrics import correct_digits, load_reference
from bench.workload import all_cases, case_id
from fatou_backend import FatouGP, find_default_gp_exe

REFERENCE_PATH = REPO / "research" / "reference" / "values.json"
DEFAULT_FORK = REPO / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"

# calibrated 2026-07-10 on the unmodified original (see --calibrate):
# required = measured worst - 5 safety digits; e|80 is a known engine outlier
# (one edge case agrees to only ~62.7 digits even in the original)
REQUIRED_DIGITS: dict[str, float] = {
    "0.8+0.4*I|80": 75.0,
    "1+I|80": 74.4,
    "10|80": 74.4,
    "2|80": 74.4,
    "2+I|80": 74.5,
    "e|80": 57.6,
    "2|200": 194.8,
    "e|200": 194.8,
}
# measured worst roundtrip on the original: 2.56e-87 (base e); frozen at *100
ROUNDTRIP_MAX = mp.mpf("3e-85")
ROUNDTRIP_BASES = ("e", "2", "10")
ROUNDTRIP_YS = ("1.5", "2.5", "3.5")
ROUNDTRIP_DPS = 80


def _groups() -> list[tuple[tuple[str, int], list]]:
    groups: dict[tuple[str, int], list] = defaultdict(list)
    for c in all_cases(deep=False):
        groups[(c.base, c.dps)].append(c)
    return sorted(groups.items(), key=lambda kv: (kv[0][1], kv[0][0]))


def measure_agreement(fork_path: Path) -> tuple[dict[str, float], list[str]]:
    """Worst agreement (digits vs frozen reference) per group; errors as report lines."""
    gp_exe = find_default_gp_exe()
    worst: dict[str, float] = {}
    errors: list[str] = []
    for (base, dps), cases in _groups():
        mp.mp.dps = dps + 40
        reference = load_reference(REFERENCE_PATH)
        group = f"{base}|{dps}"
        with FatouGP(gp_exe=gp_exe, fatou_gp=fork_path, dps=dps,
                     looplim=max(35, dps - 20), state_cache=False) as gp:
            try:
                results = gp.eval_batch(base, [f"{c.kind}({c.arg})" for c in cases])
            except RuntimeError as exc:
                errors.append(f"FAIL agree {group}: gp error: {exc}")
                continue
        worst[group] = min(correct_digits(v, reference[case_id(c)], cap=dps)
                           for c, v in zip(cases, results))
    return worst, errors


def measure_roundtrip(fork_path: Path) -> tuple[dict[str, mp.mpf], list[str]]:
    """Worst |sexp(slog(y)) - y| per base; errors as report lines."""
    gp_exe = find_default_gp_exe()
    mp.mp.dps = ROUNDTRIP_DPS + 40
    worst: dict[str, mp.mpf] = {}
    errors: list[str] = []
    for base in ROUNDTRIP_BASES:
        with FatouGP(gp_exe=gp_exe, fatou_gp=fork_path, dps=ROUNDTRIP_DPS,
                     state_cache=False) as gp:
            try:
                residuals = gp.roundtrip_residuals(base, [mp.mpf(y) for y in ROUNDTRIP_YS])
            except RuntimeError as exc:
                errors.append(f"FAIL roundtrip {base}: gp error: {exc}")
                continue
        worst[base] = max(abs(r) for r in residuals)
    return worst, errors


def run_gate(fork_path: Path) -> tuple[bool, list[str]]:
    report: list[str] = []
    ok = True

    agreement, agree_errors = measure_agreement(fork_path)
    ok = ok and not agree_errors
    report.extend(agree_errors)
    for group, required in REQUIRED_DIGITS.items():
        if group not in agreement:
            continue  # gp error already reported
        line_ok = agreement[group] >= required
        ok = ok and line_ok
        report.append(f"{'PASS' if line_ok else 'FAIL'} agree {group}: "
                      f"worst {agreement[group]:.1f} digits (required {required})")

    roundtrip, rt_errors = measure_roundtrip(fork_path)
    ok = ok and not rt_errors
    report.extend(rt_errors)
    for base, worst_res in roundtrip.items():
        line_ok = worst_res < ROUNDTRIP_MAX
        ok = ok and line_ok
        report.append(f"{'PASS' if line_ok else 'FAIL'} roundtrip {base}: "
                      f"worst {mp.nstr(worst_res, 6)} (max {mp.nstr(ROUNDTRIP_MAX, 3)})")

    return ok, report


def main() -> None:
    parser = argparse.ArgumentParser(description="fatou_fork.gp correctness gate")
    parser.add_argument("--fork", default=str(DEFAULT_FORK))
    parser.add_argument("--calibrate", action="store_true",
                        help="print measured worsts on the given fork (one-time use)")
    args = parser.parse_args()
    fork = Path(args.fork)
    if args.calibrate:
        agreement, errors = measure_agreement(fork)
        roundtrip, rt_errors = measure_roundtrip(fork)
        for line in errors + rt_errors:
            print(line)
        for group, worst in agreement.items():
            print(f"agree {group}: worst {worst:.2f} digits")
        for base, worst in roundtrip.items():
            print(f"roundtrip {base}: worst {mp.nstr(worst, 8)}")
        return
    ok, report = run_gate(fork)
    print("\n".join(report))
    print("GATE:", "PASS" if ok else "FAIL")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
