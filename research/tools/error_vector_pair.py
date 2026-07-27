"""Run an error-vector PAIR and report how many digits it proves.

The only verification this project accepts (METHODS §1): compute the same value
twice with decorrelated settings; where the two agree to d digits, the worse
run's true error is at most ~10^-d. That is how sexp_e(0.5) is proven to 972
digits, from a dps 1020/1033 pair.

Decorrelation here is by working precision, which changes the whole trajectory
-- grid sizes, the precision ladder, the number of iterations -- not just the
stopping point. That is stronger than varying nlim alone.

HOW BIG A GAP THE PAIR NEEDS (measured 2026-07-27, before committing 10 h to a
run). Agreement between the pair must not exceed either run's true accuracy:

    dps 300/304  (1.3%)  bit-IDENTICAL values -- the pair proves nothing
    dps 300/315  (5.0%)  agreement 302.2 vs truth 302.0     sound
    dps 520/526  (1.2%)  agreement 514.2 vs truth 508.7     over-claims 5.5
    dps 520/548  (5.4%)  agreement 514.2 vs truth 508.7     over-claims 5.5

Two lessons. First, a small gap can collapse entirely: at dps 300 a 1.3% gap
produced the same value twice, which is exactly the correlated-error artifact
that invalidated the v1 references (METHODS §1). Use >= 5%.

Second, and independent of the gap: at dps 520 both pairs agree with each other
5.5 digits BETTER than either agrees with truth. They share a residual error and
differ only beyond it, so the raw agreement systematically over-states accuracy.
That is why proven_digits subtracts a margin instead of reporting the agreement.

    python research/tools/error_vector_pair.py --dps 1250 1265 --pin

Deliberately does NOT write into research/reference/: that directory is
immutable (anti-gaming), and extending the proven ladder is a decision for a
human, not a side effect of a script. Results land in a candidate JSON.

Long runs are the point, so the script is built to survive them: every result
is flushed to disk the moment it exists, and a bracket run at a cheap precision
is taken before and after so that a machine which degraded mid-way is visible
afterwards rather than assumed away.
"""
from __future__ import annotations

import argparse
import json
import re as _re
import subprocess
import sys
import time
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from digits_vs_reference import pin  # verified affinity helper

REPO = Path(__file__).resolve().parents[2]
FORK = REPO / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
GP = r"C:\Program Files (x86)\Pari64-2-17-3\gp.exe"
REF = REPO / "research" / "reference" / "values.json"


def run(dps: int, parisize_gb: int) -> tuple[str, float, float]:
    """One converged sexp_e(0.5) at `dps`. Returns (value, gp seconds, wall)."""
    script = "\n".join([
        f"default(parisizemax, {parisize_gb * 1024**3});",
        f"default(realprecision, {dps});",
        f'read("{FORK.as_posix()}");',
        f"default(realprecision, {dps});",
        "quietmode=1;",
        "gettime();",
        "sexpinit(exp(1),0,0,0);",
        'print("MS ", gettime());',
        'print("VAL ", sexp(0.5));',
        "quit;",
    ]) + "\n"
    t0 = time.time()
    proc = subprocess.run([GP, "-q", "-f"], input=script, capture_output=True,
                          text=True, encoding="utf-8", errors="replace", cwd=str(REPO))
    wall = time.time() - t0
    m = _re.search(r"VAL\s+([-\d.eE ]+)", proc.stdout)
    if not m:
        raise SystemExit(f"dps {dps} produced no value:\n"
                         f"{proc.stdout[-2000:]}\n{proc.stderr[-1000:]}")
    ms = _re.search(r"MS\s+(\d+)", proc.stdout)
    value = m.group(1).strip().replace(" E", "e").replace(" ", "")
    return value, (int(ms.group(1)) / 1000 if ms else -1.0), wall


def agreement_digits(a: str, b: str) -> float:
    """Digits of agreement, capped at what was actually printed.

    Identical strings mean the two runs did not decorrelate at all, which is
    NO information about accuracy -- not infinite information. Returning a
    large number there is exactly how the v2 reference disaster read spurious
    agreement as proof, so it returns -1.0 instead and the caller must handle
    it.
    """
    printed = min(sum(c.isdigit() for c in a), sum(c.isdigit() for c in b))
    with mp.workdps(printed + 50):
        x, y = mp.mpf(a), mp.mpf(b)
        if x == y:
            return -1.0
        d = float(-mp.log10(abs(x - y) / max(abs(x), mp.mpf(1))))
    return min(d, float(printed))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dps", type=int, nargs=2, required=True,
                    help="the two working precisions, e.g. 1250 1265")
    ap.add_argument("--bracket-dps", type=int, default=520,
                    help="cheap run taken before and after, to expose machine drift")
    ap.add_argument("--parisize-gb", type=int, default=12)
    ap.add_argument("--pin", nargs="?", const="0xfff", default=None)
    ap.add_argument("--out", default="research/tools/_pair_candidate.json")
    args = ap.parse_args()

    if args.pin is not None:
        print(pin(int(args.pin, 0)), flush=True)

    out = REPO / args.out
    record: dict = {"dps_pair": args.dps, "runs": {}, "bracket": {}}

    def save() -> None:
        out.write_text(json.dumps(record, indent=1), encoding="utf-8")

    def bracket(tag: str) -> None:
        value, secs, wall = run(args.bracket_dps, args.parisize_gb)
        record["bracket"][tag] = {"dps": args.bracket_dps, "seconds": secs, "wall": wall}
        print(f"bracket {tag:6s} dps {args.bracket_dps}: {secs:8.2f} s", flush=True)
        save()

    bracket("before")
    for dps in args.dps:
        print(f"-- dps {dps} starting (this is the long one) --", flush=True)
        value, secs, wall = run(dps, args.parisize_gb)
        record["runs"][str(dps)] = {"value": value, "seconds": secs, "wall": wall,
                                    "decimals": sum(c.isdigit() for c in value)}
        print(f"   dps {dps}: {secs:9.2f} s ({secs/3600:.2f} h), "
              f"{sum(c.isdigit() for c in value)} decimals", flush=True)
        save()
    bracket("after")

    a, b = (record["runs"][str(d)]["value"] for d in args.dps)
    proven = agreement_digits(a, b)
    record["agreement_digits"] = proven
    # MEASURED safety margin. At dps 520 two decorrelated runs agree with each
    # other to 514.2 digits while each is only 508.7 digits from truth -- they
    # share a residual error and differ only beyond it. The raw agreement of a
    # pair therefore over-states accuracy by ~5.5 digits at that tier. 10 is
    # that, doubled, and floored: a proof is never rounded up.
    record["over_claim_margin"] = 10
    record["proven_digits"] = float(int(proven)) - 10 if proven > 0 else -1.0

    # cross-check against the existing ladder, which is proven to 972
    payload = json.loads(REF.read_text(encoding="utf-8"))
    ref = payload["values"]["sexp|e|0.5|1000"]["real"]
    record["vs_existing_reference_digits"] = agreement_digits(a, ref)

    lo = record["bracket"]["before"]["seconds"]
    hi = record["bracket"]["after"]["seconds"]
    record["machine_drift"] = hi / lo if lo else None
    save()

    print(f"\nagreement between the two runs : {proven:.1f} digits")
    print(f"  -> PROVEN                    : {int(proven)} digits")
    print(f"existing ladder tops out at    : 972 digits")
    print(f"agreement with it              : {record['vs_existing_reference_digits']:.1f}")
    print(f"machine drift across the pair  : {record['machine_drift']:.3f}x")
    print(f"\ncandidate written to {out}")
    print("NOT added to research/reference/ -- that is a human decision.")


if __name__ == "__main__":
    main()
