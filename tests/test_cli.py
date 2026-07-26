"""Smoke tests for the two shipped console scripts.

`fatou-cli` and `mixed-phase-cli` are both `[project.scripts]` entry points and
both are advertised in the README, and neither had a single test. That is how
`--values 0.1` shipped: the argument went through `ast.literal_eval`, became a
Python float, and the engine was asked about a number ~5e-18 away from the one
the user typed -- while the CLI printed dps-24 digits of the answer. It went
unnoticed because the documented example uses 0.5, which is binary-exact.

Kept at dps 40 so the whole file runs in a couple of seconds.
"""
from __future__ import annotations

import io
import sys
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend import cli


def run_cli(*argv: str) -> str:
    out = io.StringIO()
    with mock.patch.object(sys, "argv", ["fatou-cli", *argv]), redirect_stdout(out):
        cli.main()
    return out.getvalue()


class ParseValueTests(unittest.TestCase):
    """The parser must not round-trip decimals through a binary float."""

    def setUp(self) -> None:
        mp.mp.dps = 60

    def test_inexact_decimals_are_parsed_exactly(self) -> None:
        for text in ("0.1", "0.3", "1.7", "0.35"):
            with self.subTest(text=text):
                got = mp.re(cli._parse_value(text))
                self.assertEqual(mp.nstr(got, 30), mp.nstr(mp.mpf(text), 30))
                # the specific failure mode: agreeing with float(text) to the
                # bitter end means it went through a double
                self.assertNotEqual(got, mp.mpf(float(text)),
                                    f"{text} was routed through a Python float")

    def test_exact_binary_decimal_still_works(self) -> None:
        self.assertEqual(mp.re(cli._parse_value("0.5")), mp.mpf("0.5"))

    def test_scientific_and_complex_forms(self) -> None:
        self.assertEqual(mp.re(cli._parse_value("1e-3")), mp.mpf("0.001"))
        z = cli._parse_value("0.5+0.25j")
        self.assertEqual(mp.re(z), mp.mpf("0.5"))
        self.assertEqual(mp.im(z), mp.mpf("0.25"))

    def test_garbage_gets_a_message_not_a_traceback(self) -> None:
        with self.assertRaises(SystemExit) as ctx:
            cli._parse_value("not-a-number")
        self.assertIn("--values", str(ctx.exception))


class CliEndToEndTests(unittest.TestCase):
    def test_sexp_matches_the_known_value(self) -> None:
        out = run_cli("sexp", "--base", "e", "--values", "0.5", "--dps", "40")
        mp.mp.dps = 60
        printed = mp.mpf(out.strip().splitlines()[0])
        # The CLI prints max(30, dps-24) significant digits, i.e. 30 here, so
        # the comparison can only be as tight as that truncation allows.
        self.assertLess(abs(printed - mp.mpf("1.6463542337511945809719240315921145182")),
                        mp.mpf("1e-28"))

    def test_inexact_argument_reaches_the_engine_intact(self) -> None:
        """End-to-end guard for the --values blocker.

        sexp_e is strictly increasing, and sexp_e(0.1) differs from
        sexp_e(0.1+5e-18) far below the printed precision -- so instead of
        chasing that, compare the CLI against the library called with the
        exact string, which is the path that was always correct.
        """
        out = run_cli("sexp", "--base", "e", "--values", "0.1", "--dps", "40")
        printed = mp.mpf(out.strip().splitlines()[0])

        from fatou_backend import FatouGP
        exact = mp.re(FatouGP(dps=40).sexp("e", "0.1"))
        self.assertLess(abs(printed - exact), mp.mpf("1e-28"),
                        "CLI disagrees with the exact-string path")

    def test_missing_values_is_a_clean_error(self) -> None:
        with self.assertRaises(SystemExit):
            run_cli("sexp", "--base", "e", "--dps", "40")


class MixedPhaseTests(unittest.TestCase):
    """The other shipped entry point: exercise its per-base ladders."""

    # The ladders reduce TOWER arguments, so they expect heights where d^t is
    # already large. Below that they legitimately go complex (log of a negative
    # intermediate): inner="2" at t=0.5 returns -0.717 + pi*i. Measured, so the
    # inputs here stay in the intended regime.
    HEIGHTS = [mp.mpf("5"), mp.mpf("20"), mp.mpf("50")]

    def test_every_inner_base_branch_is_reachable(self) -> None:
        """Three hand-written ladders plus a raise, none of them exercised."""
        from fatou_backend import mixed_phase
        mp.mp.dps = 40
        for inner in ("2", "e", "10"):
            with self.subTest(inner=inner):
                got = mixed_phase.reduction_argument("e", inner, self.HEIGHTS)
                self.assertEqual(len(got), len(self.HEIGHTS))
                for v in got:
                    self.assertTrue(mp.isfinite(v),
                                    f"non-finite reduction for inner base {inner}")
                    self.assertEqual(mp.im(mp.mpc(v)), 0,
                                     f"inner base {inner} left the real axis "
                                     f"inside its intended regime")

    def test_unsupported_inner_base_raises(self) -> None:
        from fatou_backend import mixed_phase
        with self.assertRaises(ValueError):
            mixed_phase.reduction_argument("e", "7", [mp.mpf("0.5")])

    def test_reduction_is_increasing_in_the_height(self) -> None:
        """Pins direction, not just finiteness: a taller tower reduces higher."""
        from fatou_backend import mixed_phase
        mp.mp.dps = 40
        for inner in ("2", "e", "10"):
            with self.subTest(inner=inner):
                got = mixed_phase.reduction_argument("e", inner, self.HEIGHTS)
                self.assertEqual(list(got), sorted(got),
                                 f"inner base {inner}: reduction not monotone "
                                 f"in the tower height")


if __name__ == "__main__":
    unittest.main()
