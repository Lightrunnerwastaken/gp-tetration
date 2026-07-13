"""Re-verify any entry of research/reference/values.json from scratch.

Every reference entry documents its provenance (engine, dps pair, error
vector) in the meta block of values.json. This script re-runs the engine
and reports the agreement, so anyone with the public repo + PARI/GP can
reproduce the digit claims:

    python research/tools/verify_reference.py --list
    python research/tools/verify_reference.py --key "sexp|e|0.5|80"
    python research/tools/verify_reference.py --key "sexp|e|0.5|500" --dps 520
    python research/tools/verify_reference.py --key "sexp|e|0.5|500" --dps 520 --pair 533

Modes:
  default   run once, compare against the stored reference value
  --pair D2 additionally run at dps D2 and report run-vs-run agreement —
            a fresh error-vector proof, independent of the stored value

Approximate runtimes with the fork (base e): dps 100 < 1 min, dps 300
~3 min, dps 520 ~25 min, dps 1020 ~6.7 h. Base 2 is ~2x slower per digit.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

import mpmath as mp

from bench.metrics import correct_digits, load_reference
from fatou_backend import FatouGP

VALUES = REPO / "research" / "reference" / "values.json"
FORK = REPO / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
ORIGINAL = REPO / "src" / "fatou_backend" / "vendor" / "fatou.gp"


def parse_key(key: str) -> tuple[str, str, str]:
    op, base, x, _tier = key.split("|")
    if op not in ("sexp", "slog"):
        raise SystemExit(f"unsupported op {op!r} in key {key!r}")
    if base.lower() == "e":
        base = "exp(1)"
    return op, base, x


def provenance(payload: dict, key: str) -> list[str]:
    notes = []
    for section, block in payload.get("meta", {}).items():
        if isinstance(block, dict) and key in block:
            entry = block[key]
            if isinstance(entry, dict) and "method" in entry:
                notes.append(f"[{section}] {entry['method']}")
    return notes


def run_engine(op: str, base: str, x: str, dps: int, args) -> mp.mpc:
    engine = ORIGINAL if args.original else FORK
    # knob conventions match the recorded provenance: the fork converges
    # via looplim = digits target (as in gate.py); the original via a
    # high iteration cap (nlim) with looplim=0.
    if args.original:
        nlim = args.nlim if args.nlim is not None else 400
        looplim = args.looplim if args.looplim is not None else 0
    else:
        nlim = args.nlim if args.nlim is not None else 30
        looplim = args.looplim if args.looplim is not None else max(35, dps - 20)
    gp = FatouGP(fatou_gp=engine, dps=dps, nlim=nlim, looplim=looplim,
                 state_cache=False)
    t0 = time.perf_counter()
    value = getattr(gp, op)(base, x)
    print(f"  {op}({base}, {x}) @ dps {dps} "
          f"[{engine.name}, nlim={nlim}, looplim={looplim}]: "
          f"{time.perf_counter() - t0:.1f} s")
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--key", help='reference key, e.g. "sexp|e|0.5|500"')
    parser.add_argument("--list", action="store_true", help="list keys and exit")
    parser.add_argument("--dps", type=int, default=None,
                        help="working precision (default: tier + 20)")
    parser.add_argument("--pair", type=int, default=None, metavar="DPS2",
                        help="second run at this dps for a fresh error-vector proof")
    parser.add_argument("--nlim", type=int, default=None,
                        help="iteration cap (default: engine-appropriate)")
    parser.add_argument("--looplim", type=int, default=None,
                        help="convergence target (default: max(35, dps-20) "
                             "for the fork, 0 for the original)")
    parser.add_argument("--original", action="store_true",
                        help="use the unmodified fatou.gp instead of the fork")
    args = parser.parse_args()

    payload = json.loads(VALUES.read_text(encoding="utf-8"))
    if args.list or not args.key:
        for key in payload["values"]:
            print(key)
        return

    if args.key not in payload["values"]:
        raise SystemExit(f"unknown key {args.key!r} (use --list)")
    op, base, x = parse_key(args.key)
    tier = int(args.key.split("|")[-1])
    dps = args.dps or tier + 20
    # set precision BEFORE parsing the reference strings, or they truncate
    mp.mp.dps = max(mp.mp.dps, dps + 20, tier + 20)
    references = load_reference(VALUES)

    for note in provenance(payload, args.key):
        print(f"stored provenance: {note}")

    value = run_engine(op, base, x, dps, args)
    agree = correct_digits(value, references[args.key], cap=dps)
    print(f"agreement with stored reference: {agree:.1f} digits")

    if args.pair:
        value2 = run_engine(op, base, x, args.pair, args)
        pair_agree = correct_digits(value, value2, cap=min(dps, args.pair))
        print(f"run-vs-run agreement (error vector): {pair_agree:.1f} digits")
        print("=> the worse run's true error is bounded by this agreement")


if __name__ == "__main__":
    main()
