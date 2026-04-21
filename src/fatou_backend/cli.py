from __future__ import annotations

import argparse
import ast

import mpmath as mp

from .gp_backend import FatouGP


def _parse_value(text: str) -> mp.mpc:
    try:
        parsed = ast.literal_eval(text)
    except Exception:
        parsed = text
    return mp.mpc(parsed)


def _print_values(values: list[mp.mpc]) -> None:
    for value in values:
        if abs(mp.im(value)) <= mp.mpf("1e-30"):
            print(mp.nstr(mp.re(value), 30))
        else:
            print(f"{mp.nstr(mp.re(value), 30)} + {mp.nstr(mp.im(value), 30)}*I")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="CLI wrapper for fatou.gp")
    parser.add_argument("command", choices=["sexp", "slog", "roundtrip", "eval"])
    parser.add_argument("--base", required=True, help="GP base expression, e.g. e, 2, 10, 1+I")
    parser.add_argument("--values", nargs="+", help="Values to evaluate")
    parser.add_argument("--expressions", nargs="+", help="Raw GP expressions for eval mode")
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--nlim", type=int, default=30)
    parser.add_argument("--nskip", type=int, default=4)
    parser.add_argument("--looplim", type=int, default=35)
    parser.add_argument("--gp-exe")
    parser.add_argument("--fatou-gp")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    mp.mp.dps = max(args.dps + 20, 80)
    gp = FatouGP(
        gp_exe=args.gp_exe,
        fatou_gp=args.fatou_gp,
        dps=args.dps,
        nlim=args.nlim,
        nskip=args.nskip,
        looplim=args.looplim,
    )
    session = gp.session(args.base)
    if args.command == "eval":
        if not args.expressions:
            raise SystemExit("--expressions is required for eval.")
        result = session.eval_batch(args.expressions)
    else:
        if not args.values:
            raise SystemExit("--values is required for sexp, slog, and roundtrip.")
        values = [_parse_value(text) for text in args.values]
        if args.command == "sexp":
            result = session.sexp_batch(values)
        elif args.command == "slog":
            result = session.slog_batch(values)
        else:
            result = session.roundtrip_residuals(values)
    _print_values(result)


if __name__ == "__main__":
    main()
