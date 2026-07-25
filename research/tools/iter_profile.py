"""exp-040 (H3): where does sexpinit's time actually go, per iteration?

The evidence H3 rests on (exp-009, 2026-07-10: "letzte 10 Iterationen je ~10.4 s")
predates ~20 keeps (FFT extraction exp-011, 256-grids exp-020, the isuperf cache
exp-024, theta quantization exp-026/027, chirp caches exp-030/031, ...). The cost
distribution must be re-measured before any mutation is designed.

Instrument v2. v1 used nlim differencing and failed (Run 20260725-010259,
INCONCLUSIVE): `re` is local to loop() and unreadable after return, and the run
saturates because loop()'s stall-exit fires before nlim binds. v2 instead reads
loop()'s own per-iteration line, which quietmode=0 already prints:

    n=loopcnt RE decimal digits, CT ctsamples, TH thsamples

No fork mutation. Each line is timestamped on arrival; if GP block-buffers when
piped, the timestamps collapse and the tool says so instead of reporting fiction
— the iteration/ctsamples curve is still valid either way.

  python research/tools/iter_profile.py --base "exp(1)" --dps 100
"""
from __future__ import annotations

import argparse
import re as _re
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
FORK = REPO / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DEFAULT_GP = r"C:\Program Files (x86)\Pari64-2-17-3\gp.exe"

_ITER = _re.compile(
    r"^(\d+)=loopcnt\s+([-\d.eE]+)\s+decimal digits,\s*(\d+)\s+ctsamples,\s*(\d+)\s+thsamples")


def _script(base: str, dps: int) -> str:
    return "\n".join([
        "default(parisizemax, 2147483648);",
        f"default(realprecision, {dps});",
        f'read("{FORK.as_posix()}");',
        "quietmode=0;",
        "gettime();",
        f"sexpinit({base},0,0,0);",
        'print("PROBE-MS ", gettime());',
        "quit;",
    ]) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--base", default="exp(1)")
    ap.add_argument("--dps", type=int, default=100)
    ap.add_argument("--gp", default=DEFAULT_GP)
    args = ap.parse_args()

    print(f"iter-profile v2: base={args.base} dps={args.dps}", flush=True)
    proc = subprocess.Popen(
        [args.gp, "-q", "-f"], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace",
        bufsize=1, cwd=str(REPO),
    )
    proc.stdin.write(_script(args.base, args.dps))
    proc.stdin.close()

    t0 = time.monotonic()
    prev_t, prev_re = t0, None
    iters, total_ms = [], None
    for line in proc.stdout:
        line = line.strip()
        m = _ITER.match(line)
        if m:
            now = time.monotonic()
            n, re_d, ct, th = int(m.group(1)), float(m.group(2)), int(m.group(3)), int(m.group(4))
            dt_ms = (now - prev_t) * 1000
            gained = (re_d - prev_re) if prev_re is not None else re_d
            iters.append((n, re_d, ct, th, dt_ms, gained))
            print(f'KERN-METRIC {{"iter": {n}, "re": {re_d:.2f}, "ctsamples": {ct}, '
                  f'"thsamples": {th}, "dt_ms": {dt_ms:.0f}, "digits_gained": {gained:.2f}}}',
                  flush=True)
            prev_t, prev_re = now, re_d
        elif line.startswith("PROBE-MS"):
            total_ms = int(line.split()[1])
    proc.wait(timeout=60)

    if not iters:
        print("RESULT NO-DATA: loop() printed no iteration lines (quietmode?)", flush=True)
        return 1

    wall_ms = (time.monotonic() - t0) * 1000
    measured = sum(r[4] for r in iters)
    buffered = measured < 0.25 * wall_ms          # arrivals collapsed -> timing unusable

    print(f"\n  iter    re   ctsamples   dt_ms  digits", flush=True)
    for n, re_d, ct, th, dt_ms, gained in iters:
        print(f"  {n:4d} {re_d:7.1f} {ct:9d} {dt_ms:8.0f} {gained:7.2f}", flush=True)

    last = iters[-1]
    print(f'\nKERN-METRIC {{"iterations": {last[0]}, "final_re": {last[1]:.2f}, '
          f'"final_ctsamples": {last[2]}, "total_ms": {total_ms}, '
          f'"buffered": {str(buffered).lower()}}}', flush=True)

    if buffered:
        print("RESULT TIMING-UNUSABLE: GP block-buffered its output (arrival deltas "
              f"sum to {measured:.0f} ms of {wall_ms:.0f} ms wall). The iteration and "
              "ctsamples curve above is still valid; per-iteration timing is not.",
              flush=True)
        return 0

    half = len(iters) // 2
    early = sum(r[4] for r in iters[:half]) / max(half, 1)
    late = sum(r[4] for r in iters[half:]) / max(len(iters) - half, 1)
    ratio = late / early if early else 0
    print(f'KERN-METRIC {{"late_vs_early_per_iter": {ratio:.2f}}}', flush=True)
    print(f"RESULT {'LATE-DOMINATED' if ratio >= 2.0 else 'FLAT'}: second half costs "
          f"{ratio:.2f}x the first per iteration (H3 expects >=2x)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
