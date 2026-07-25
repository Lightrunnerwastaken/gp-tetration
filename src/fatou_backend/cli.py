from __future__ import annotations

import argparse
import ast

import mpmath as mp

from .gp_backend import FatouGP


def _parse_value(text: str) -> mp.mpc:
    """Parse a --values argument WITHOUT going through a Python float.

    The old form ran ast.literal_eval() first, so any decimal that is not
    exactly representable in binary -- 0.1, 0.3, 1.7, i.e. most things a user
    types -- became a float. The engine was then asked about a different
    argument than the user wrote, off by ~5e-18, and the answer was printed to
    dps-24 digits with no hint that only ~17 of them meant anything. It went
    unnoticed because the documented example uses 0.5, which is binary-exact.

    mpmath parses decimal strings exactly at the current precision, so try that
    first and keep literal_eval only for Python-syntax complex forms.
    """
    try:
        return mp.mpc(mp.mpf(text))
    except (ValueError, TypeError):
        pass
    try:
        return mp.mpc(text)
    except (ValueError, TypeError):
        pass
    try:
        return mp.mpc(ast.literal_eval(text))
    except Exception as exc:
        raise SystemExit(
            f"--values: cannot parse {text!r} as a number. Use a decimal "
            f"(0.5, 1e-3), or Python complex syntax (0.5+0.25j)."
        ) from exc


def _print_values(values: list[mp.mpc], digits: int = 30) -> None:
    for value in values:
        if abs(mp.im(value)) <= mp.mpf("1e-30"):
            print(mp.nstr(mp.re(value), digits))
        else:
            print(f"{mp.nstr(mp.re(value), digits)} + {mp.nstr(mp.im(value), digits)}*I")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="CLI wrapper for fatou.gp")
    parser.add_argument("command", choices=["sexp", "slog", "roundtrip", "eval"])
    parser.add_argument("--base", required=True, help="GP base expression, e.g. e, 2, 10, 1+I")
    parser.add_argument("--values", nargs="+", help="Values to evaluate")
    parser.add_argument("--expressions", nargs="+", help="Raw GP expressions for eval mode")
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--nlim", type=int, default=None,
                        help="iteration cap (default: scales with dps)")
    parser.add_argument("--nskip", type=int, default=4)
    parser.add_argument("--looplim", type=int, default=None,
                        help="convergence target (default: full precision)")
    parser.add_argument("--gp-exe")
    parser.add_argument("--fatou-gp",
                        help='engine file, or the shortcuts "fork" / "original"')
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
    # Deliberately conservative: the fork delivers >=0.95*dps true digits, so
    # dps-24 under-prints above ~dps 500. Kept as-is because the 0.95 rule is
    # measured only at dps 300/400/520 and this display floor must never
    # over-print at tiers nobody has checked.
    _print_values(result, digits=max(30, args.dps - 24))


if __name__ == "__main__":
    main()
