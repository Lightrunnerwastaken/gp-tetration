"""Measure how many TRUE digits an engine delivers, against the proven references.

Comparing two engines against each other only shows where they diverge, not
which one is right. research/reference/values.json holds error-vector-proven
values (sexp_e(0.5) to 972 digits), so the honest question -- "did this
mutation cost real digits?" -- is answered here, not by the gate (whose e-group
thresholds are a low legacy of the v2 correction) and not by mutual agreement.

    python research/tools/digits_vs_reference.py --dps 300 \
        --engine src/fatou_backend/vendor/fatou_fork.gp \
        --engine research/tools/_exp044_1_60.gp
"""
from __future__ import annotations

import argparse
import ctypes
import json
import re as _re
import subprocess
import sys
import time
from pathlib import Path

import mpmath as mp


def pin(mask: int) -> str:
    """Pin this process (and every gp.exe it spawns) to a fixed core set.

    On a hybrid-core CPU the Windows scheduler will migrate a long-running
    background job onto the E-cores as soon as foreground apps appear. Measured
    on an i7-12700H, dps 300, same engine, digit-identical output: 51 s pinned
    to the P-cores vs 85 s unpinned -- a 1.68x swing against a calibrated noise
    threshold of 6.8%. Any timing comparison across precisions is meaningless
    without this. Children inherit the mask, so gp.exe is covered.

    0xfff = logical processors 0-11 = the six P-cores of a 12700H. Adjust for
    other topologies; `wmic cpu get NumberOfCores,NumberOfLogicalProcessors`.
    """
    if not sys.platform.startswith("win"):
        raise SystemExit("--pin is Windows-only; drop it or add a POSIX path")
    k32 = ctypes.windll.kernel32
    # restype MUST be set. Without it ctypes hands back the pseudo-handle as a
    # signed 32-bit int, which widens to 0x00000000FFFFFFFF instead of
    # 0xFFFFFFFFFFFFFFFF -- an invalid handle, so the call fails and the
    # process silently keeps running on every core. That exact bug produced a
    # 20-minute "pinned" measurement that was not pinned at all.
    k32.GetCurrentProcess.restype = ctypes.c_void_p
    k32.SetProcessAffinityMask.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
    k32.GetProcessAffinityMask.argtypes = [ctypes.c_void_p,
                                           ctypes.POINTER(ctypes.c_size_t),
                                           ctypes.POINTER(ctypes.c_size_t)]
    h = k32.GetCurrentProcess()
    if not k32.SetProcessAffinityMask(h, mask):
        raise SystemExit(f"--pin: SetProcessAffinityMask({mask:#x}) failed, "
                         f"err {ctypes.get_last_error()}")
    # Read it back. A pin that quietly does nothing is worse than no pin: it
    # makes an untrustworthy measurement look trustworthy.
    got, sysmask = ctypes.c_size_t(), ctypes.c_size_t()
    if not k32.GetProcessAffinityMask(h, ctypes.byref(got), ctypes.byref(sysmask)):
        raise SystemExit("--pin: could not read the affinity back")
    if got.value != mask:
        raise SystemExit(f"--pin: asked for {mask:#x}, got {got.value:#x} "
                         f"(system mask {sysmask.value:#x}) -- refusing to "
                         f"report timings that are not actually pinned")
    return f"pinned to mask {mask:#x} (verified)"

REPO = Path(__file__).resolve().parents[2]
GP = r"C:\Program Files (x86)\Pari64-2-17-3\gp.exe"
REF = REPO / "research" / "reference" / "values.json"


def reference(key: str) -> tuple[mp.mpf, float | None]:
    """Return the reference value AND the number of digits actually PROVEN.

    values.json stores 1000 decimals under sexp|e|0.5|1000, but only 972 of
    them are error-vector proven -- the rest are the producing engine's own
    unverified tail. Agreement past the proven ceiling is partly self-
    agreement and must not be reported as true digits. (This tool printed
    "992.0 true digits" for a dps-1020 run before the ceiling was enforced.)
    """
    payload = json.loads(REF.read_text(encoding="utf-8"))
    e = payload["values"][key]
    proven: float | None = None

    def scan(o) -> None:
        nonlocal proven
        if isinstance(o, dict):
            for k, v in o.items():
                if k == "proven_digits" and isinstance(v, dict) and key in v:
                    proven = float(v[key])
                else:
                    scan(v)
        elif isinstance(o, list):
            for v in o:
                scan(v)

    scan(payload.get("meta", {}))
    return mp.mpf(e["real"]), proven


def run(engine: Path, dps: int, base: str, arg: str, extra: list[str]) -> tuple[str, int]:
    script = "\n".join([
        "default(parisizemax, 8589934592);",
        f"default(realprecision, {dps});",
        f'read("{engine.as_posix()}");',
        f"default(realprecision, {dps});",
        "quietmode=1;",
        *[f"{e};" for e in extra],
        f"rbase = {base};",
        "gettime();",
        "sexpinit(rbase,0,0,0);",
        'print("RMS ", gettime());',
        f'print("RVAL ", sexp({arg}));',
        "quit;",
    ]) + "\n"
    p = subprocess.run([GP, "-q", "-f"], input=script, capture_output=True,
                       text=True, encoding="utf-8", errors="replace", cwd=str(REPO))
    v = _re.search(r"RVAL\s+([-\d.eE ]+)", p.stdout)
    m = _re.search(r"RMS\s+(\d+)", p.stdout)
    if not v:
        raise SystemExit(p.stdout[-1200:] + p.stderr[-500:])
    return v.group(1).strip().replace(" E", "e").replace(" ", ""), int(m.group(1)) if m else -1


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", action="append", required=True)
    ap.add_argument("--dps", type=int, nargs="+", default=[300])
    ap.add_argument("--base", default="exp(1)")
    ap.add_argument("--arg", default="0.5")
    ap.add_argument("--key", default="sexp|e|0.5|1000")
    ap.add_argument("--set", action="append", default=[])
    ap.add_argument("--pin", nargs="?", const="0xfff", default=None,
                    help="pin to a core mask (default 0xfff = P-cores of a "
                         "12700H). Required for any cross-precision comparison.")
    ap.add_argument("--repeat", type=int, default=1,
                    help="run each engine N times and report the best (the "
                         "best run is the one least disturbed by other load)")
    args = ap.parse_args()

    if args.pin is not None:
        print(pin(int(args.pin, 0)))

    mp.mp.dps = 1100
    ref, proven = reference(args.key)
    print(f"reference {args.key} loaded ({mp.mp.dps} dps working"
          + (f", proven to {proven:.0f} digits)" if proven else ", proven depth unknown)"))
    for dps in args.dps:
        print(f"\n-- dps {dps} --")
        for e in args.engine:
            eng = Path(e) if Path(e).is_absolute() else REPO / e
            t0 = time.time()
            best: tuple[str, int] | None = None
            for _ in range(max(1, args.repeat)):
                cand = run(eng, dps, args.base, args.arg, args.set)
                if best is None or cand[1] < best[1]:
                    best = cand
            val, ms = best  # type: ignore[misc]
            got = mp.mpf(val)
            err = abs(got - ref)
            digits = float(-mp.log10(err / max(abs(ref), mp.mpf(1)))) if err else float(dps)
            if proven is not None and digits > proven:
                # Saturated the reference: everything past `proven` is the
                # producing engine's unverified tail, so report the ceiling.
                shown = (f">={proven:8.1f}   (reference exhausted; "
                         f"raw agreement {digits:.1f})")
            else:
                shown = f"{digits:8.1f}"
            print(f"  {eng.name:26s} {ms/1000:8.2f}s   TRUE digits {shown}"
                  f"   (wall {time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
