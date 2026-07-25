"""exp-039 pilot: does a converged lower-dps sexpinit state warm-start a higher-dps init?

Static reading of the fork says no: sexpinit -> loop -> initsch recomputes every
state variable at the CURRENT realprecision (precis, L=fixedk, sfunczero, z0h/z0l,
the Schroeder series) and then loop() resets ct/re/n. A loaded state would be
overwritten. This probe checks that empirically in the most favourable possible
setting: the low-dps state is still live IN THE SAME PROCESS (no writebin
round-trip, no serialization loss). If even that does not speed up the high-dps
init, no state-dump scheme can.

  warm : realprecision=DPS_LO, sexpinit  ->  realprecision=DPS_HI, sexpinit (timed)
  cold : realprecision=DPS_HI, sexpinit (timed)          [fresh process]

Emits KERN-METRIC lines for the Kern runner. No fork mutation (active-plan rule).
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
FORK = REPO / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DEFAULT_GP = r"C:\Program Files (x86)\Pari64-2-17-3\gp.exe"

_MS = re.compile(r"^PROBE-MS\s+(\d+)\s*$", re.MULTILINE)


def _script(base: str, dps_lo: int, dps_hi: int, warm: bool) -> str:
    fork = FORK.as_posix()
    lines = [
        "default(parisizemax, 2147483648);",
        'quietmode=1;',
    ]
    if warm:
        lines += [
            f"default(realprecision, {dps_lo});",
            f'read("{fork}");',
            "quietmode=1;",
            f"sexpinit({base},0,0,0);",
            # low-dps state now live in this process
            f"default(realprecision, {dps_hi});",
            "gettime();",
            f"sexpinit({base},0,0,0);",
        ]
    else:
        lines += [
            f"default(realprecision, {dps_hi});",
            f'read("{fork}");',
            "quietmode=1;",
            "gettime();",
            f"sexpinit({base},0,0,0);",
        ]
    lines += [
        'print("PROBE-MS ", gettime());',
        'print("PROBE-ANCHOR ", sexp(0.5));',
        "quit;",
    ]
    return "\n".join(lines) + "\n"


def _run(gp_exe: str, script: str, timeout: float) -> tuple[int, str]:
    proc = subprocess.run(
        [gp_exe, "-q", "-f"], input=script, cwd=str(REPO),
        capture_output=True, text=True, encoding="utf-8", errors="replace",
        timeout=timeout,
    )
    out = proc.stdout + proc.stderr
    match = _MS.search(out)
    if match is None:
        raise RuntimeError(f"probe produced no timing marker:\n{out[-1500:]}")
    return int(match.group(1)), out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base", default="2")
    ap.add_argument("--dps-lo", type=int, default=60)
    ap.add_argument("--dps-hi", type=int, default=100)
    ap.add_argument("--gp", default=DEFAULT_GP)
    ap.add_argument("--timeout", type=float, default=1800.0)
    args = ap.parse_args()

    print(f"probe: base={args.base} {args.dps_lo} -> {args.dps_hi} dps", flush=True)

    cold_ms, cold_out = _run(args.gp, _script(args.base, args.dps_lo, args.dps_hi, False),
                             args.timeout)
    print(f'KERN-METRIC {{"phase": "cold", "init_ms": {cold_ms}}}', flush=True)

    warm_ms, warm_out = _run(args.gp, _script(args.base, args.dps_lo, args.dps_hi, True),
                             args.timeout)
    print(f'KERN-METRIC {{"phase": "warm", "init_ms": {warm_ms}}}', flush=True)

    speedup = (cold_ms - warm_ms) / cold_ms if cold_ms else 0.0
    print(f'KERN-METRIC {{"cold_ms": {cold_ms}, "warm_ms": {warm_ms}, '
          f'"saving_frac": {speedup:.4f}}}', flush=True)

    # H2 asks for >= 30% fewer init iterations; time is the honest proxy here
    verdict = "SUPPORTED" if speedup >= 0.30 else "REFUTED"
    print(f"RESULT {verdict}: cold {cold_ms} ms vs warm {warm_ms} ms "
          f"({speedup*100:.1f}% saved; H2 needs >=30%)", flush=True)

    for label, out in (("cold", cold_out), ("warm", warm_out)):
        anchor = [l for l in out.splitlines() if l.startswith("PROBE-ANCHOR")]
        if anchor:
            print(f"  {label} anchor: {anchor[0][13:60]}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
