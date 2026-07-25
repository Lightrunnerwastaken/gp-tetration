"""Run the sfunc branch probe at a given dps/base and report the split.

    python research/tools/branchprobe.py --dps 200 --base "exp(1)"

Reports: direct vs theta samples, full ct-Horner count, and the two
decision margins (imin = min |imag(y2)|, tmin = min relative period-tie
margin). Small margins would forbid taking those decisions at reduced
precision; O(1) margins allow it.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
GP_EXE = r"C:\Program Files (x86)\Pari64-2-17-3\gp.exe"


def run(dps: int, base: str, parisize: int = 4_000_000_000, quiet: bool = True) -> dict:
    driver = (HERE / "branchprobe_run.gp").read_text(encoding="utf-8")
    head = f'brdps={dps}; brbasestr="{base}"; brquiet={1 if quiet else 0};\n'
    tmp = HERE / f"_bp_{dps}.gp"
    tmp.write_bytes((head + driver).encode("utf-8"))
    t0 = time.time()
    proc = subprocess.run(
        [GP_EXE, "-q", "--default", f"parisize={parisize}", str(tmp)],
        capture_output=True, text=True, cwd=str(HERE),
    )
    wall = time.time() - t0
    out = proc.stdout
    res: dict = {"dps": dps, "base": base, "wall_s": round(wall, 1)}
    m = re.search(r"BRPROBE dps=(\d+) time_ms=(\d+)", out)
    if m:
        res["init_ms"] = int(m.group(2))
    m = re.search(r"direct=(\d+) up=(\d+) dn=(\d+) horner=(\d+) shift=(\d+)", out)
    if m:
        res |= {
            "direct": int(m.group(1)), "up": int(m.group(2)), "dn": int(m.group(3)),
            "horner": int(m.group(4)), "shift": int(m.group(5)),
        }
    m = re.search(r"imin=([-\d.eE]+) tmin=([-\d.eE]+) thetacalls=(\d+)", out)
    if m:
        res |= {"imin": float(m.group(1)), "tmin": float(m.group(2)),
                "theta": int(m.group(3))}
    m = re.search(r"BRPROBE anchor=([-\d.eE]+)", out)
    if m:
        res["anchor"] = m.group(1)[:40]
    if "direct" not in res:
        res["raw_tail"] = out[-600:]
        res["stderr"] = proc.stderr[-400:]
    tmp.unlink(missing_ok=True)
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dps", type=int, nargs="+", default=[200])
    ap.add_argument("--base", default="exp(1)")
    ap.add_argument("--parisize", type=int, default=4_000_000_000)
    args = ap.parse_args()
    for dps in args.dps:
        r = run(dps, args.base, args.parisize)
        tot = r.get("direct", 0) + r.get("up", 0) + r.get("dn", 0)
        if tot:
            r["theta_pct"] = round(100 * (r["up"] + r["dn"]) / tot, 1)
        print(r)


if __name__ == "__main__":
    main()
