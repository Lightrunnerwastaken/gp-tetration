"""Hexation in the fork: regular below eta5 (exp-077), Kneser above (exp-079).

For real 1 < b < eta5 = 1.7373559... the pentation crosses y = x twice and
the lower crossing a is attracting, 0 < P'(a) < 1. There the hexation is the
regular iteration -- the same construction the fork uses for tetration below
eta: hex(z) = Psi_H(kappa * mu^z), hex(0) = 1, hex(z+1) = pent(hex(z)),
hex(z) -> a. It carries the engine's precision.

Above eta5 (2, e, 10, ...) the pentation has no real fixed point. hexinit
then builds the Kneser hexation at the complex fixed point p of the pentation
(exp-079): one least-squares problem for a real Taylor series around a
planned centre on the axis (exp-081) plus the theta coefficients, solved on
two arcs so that the two solutions bound the error (exp-084..087). Sheldon Levenstein's matrix approximation (~8
digits) is only the fallback now and is checked against the Kneser value.

Cost: ~35 s for the regular part with a warm state cache, ~30 s for base 2
above eta5; base e only with FATOU_BACKEND_RUN_SLOW=1.
"""
from __future__ import annotations

import os
import sys
import unittest
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend import FatouGP, find_default_gp_exe

DPS = 40
DEEP = 60

# hexinit(B, n, skip): skip=1 reuses the engine the worker already initialised.
INIT = "hexinit(sexp(1),0,1)"


def _hex(engine: FatouGP, base: str, exprs: list[str]) -> list[mp.mpc]:
    """hexinit in the same worker, then the expressions; drops its estimate."""
    return engine.eval_batch(base, [INIT] + exprs)[1:]


class HexationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        mp.mp.dps = 120
        cls.fork = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=DPS)
        cls.deep = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=DEEP)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.fork.close()
        cls.deep.close()

    # Measured at dps 40 with fresh engines: worst residual 1.5e-47 (b = 1.73),
    # 2.5e-48 (b = 1.7), 1e-57 (b = 1.2). Threshold: the engine's ~dps digits.
    TOL = mp.mpf(f"1e-{DPS}")

    def _assert_identities(self, base: str) -> None:
        exprs = ["hexreg", "hex(0)-1", "hex(1)-sexp(1)", "hex(2)-pent(sexp(1))",
                 "hex(-1)", "hex(-2)+1", "hex(-3)+2", "ihex(hex(1/3))-1/3",
                 "pent(hexa)-hexa"]
        vals = _hex(self.fork, base, exprs)
        self.assertEqual(vals[0], 1, f"base {base} should take the regular path")
        for expr, val in zip(exprs[1:], vals[1:]):
            with self.subTest(base=base, check=expr):
                self.assertLess(abs(val), self.TOL,
                                f"base {base}: {expr} = {mp.nstr(val, 5)}")

    def test_identities_far_below_eta(self) -> None:
        self._assert_identities("1.2")

    def test_identities_between_eta_and_eta4(self) -> None:
        self._assert_identities("1.5")

    def test_identities_close_to_eta5(self) -> None:
        # P'(a) = 0.449; at 1.73 it is 0.679 and still passes, measured.
        self._assert_identities("1.7")

    def test_hexation_tends_to_the_attracting_fixed_point(self) -> None:
        mu, a, h40 = _hex(self.fork, "1.7", ["hexmu", "hexa", "hex(40)"])
        self.assertTrue(0 < mp.re(mu) < 1, f"P'(a) = {mu}")
        # hex(40) - a ~ kappa * mu^40 ~ 1e-14: convergence, not error
        self.assertLess(abs(h40 - a), mp.mpf("1e-12"))
        self.assertAlmostEqual(float(mp.re(a)), 2.3876372817, places=9)

    # Error vector, threshold at the engine's documented floor -- see
    # test_pentation.py FLOOR for why a warm state cache costs digits.
    FLOOR = mp.mpf(10) ** (-0.95 * DPS)

    def test_error_vector_against_deeper_run(self) -> None:
        # Measured with fresh engines: 1e-54..1e-57 (1.2), 1e-46..1e-51 (1.7).
        exprs = ["hex(1/2)", "hex(-1/2)", "hex(3/2)", "hex(-5/2)", "hexa", "ihex(1.1)"]
        for base in ("1.2", "1.7"):
            shallow = _hex(self.fork, base, exprs)
            deep = _hex(self.deep, base, exprs)
            for expr, a, b in zip(exprs, shallow, deep):
                with self.subTest(base=base, value=expr):
                    rel = abs(a - b) / abs(b)
                    self.assertLess(rel, self.FLOOR,
                                    f"base {base}: {expr} moves by {mp.nstr(rel, 5)} "
                                    f"between dps {DPS} and {DEEP}")



# -- above eta5: Kneser at the complex fixed point of the pentation (exp-079) --
# Run at dps 20 (working precision 38): base 2 takes ~14 s there and ~50 s at
# dps 45 (two least-squares solves each).
KN_DPS = 20

# Base 2, hex(3/2): dps 45 run, error vector 43.4 digits; the dps-20 run
# (32.6 digits) agrees to 2e-32. Sheldon's matrix_h gives 2.6728939660265...,
# 7e-9 off -- matching the 8.1 digits it reports.
KN_REF_2 = "2.6728939593426194891907964716798365760"
# Base e, hex(3/2): dps 45 run, error vector 42.2 digits (exp-081; before
# that base e stopped near 26 digits -- FE images outside the arc). The dps-20
# run (30.8 digits) agrees to 1e-30.
KN_REF_E = "9.6455039399646724496319602084019208883"
# Base 10, hex(1/2): dps 45, error vector 43.6 digits; dps 20 agrees to 2e-33.
KN_REF_10_HALF = "2.0378707129541331949023918286127848393"

RUN_SLOW = os.getenv("FATOU_BACKEND_RUN_SLOW") == "1"


class KneserHexationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        mp.mp.dps = 120
        cls.fork = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=KN_DPS)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.fork.close()

    TOL = mp.mpf("1e-25")

    def _check(self, base: str, ref: str, with_p2: bool, ref_tol: str = "1e-25") -> None:
        exprs = ["hexreg", "hexerr", "hex(3/2)", "hex(0)-1", "hex(1)-sexp(1)",
                 "hex(-1)", "hex(-2)+1", "ihex(hex(1/3))-1/3",
                 "ihex(hex(3/10+I/5))-(3/10+I/5)"]
        if with_p2:
            exprs.append("hex(2)-pent(sexp(1))")
        vals = _hex(self.fork, base, exprs)
        self.assertEqual(vals[0], 2, f"base {base} should take the Kneser path")
        self.assertGreater(mp.re(vals[1]), 25, f"base {base}: error vector says {vals[1]} digits")
        rel = abs(mp.re(vals[2]) - mp.mpf(ref)) / mp.mpf(ref)
        self.assertLess(rel, mp.mpf(ref_tol), f"base {base}: hex(3/2) is {mp.nstr(rel, 5)} off the reference")
        for expr, val in zip(exprs[3:], vals[3:]):
            with self.subTest(base=base, check=expr):
                self.assertLess(abs(val), self.TOL, f"base {base}: {expr} = {mp.nstr(val, 5)}")

    def test_base_2(self) -> None:
        # P(2) = 2^^^2 = 4, so hex(2) - pent(2) is a real check here
        self._check("2", KN_REF_2, with_p2=True)

    def test_sheldon_agrees_to_his_eight_digits(self) -> None:
        """Levenstein, Tetration Forum tid 1263 (2020): 2^^^^1.5 ~ 2.6729."""
        vals = self.fork.eval_batch("2", ["pentinit(sexp(1),1)", "hexreg=0", "matrix_h(256)", "hex(3/2)"])
        diff = abs(mp.re(vals[3]) - mp.mpf(KN_REF_2))
        self.assertLess(diff, mp.mpf("1e-7"))
        self.assertGreater(diff, mp.mpf("1e-12"), "Sheldon's path is suddenly exact -- which path ran?")

    @unittest.skipUnless(RUN_SLOW, "Set FATOU_BACKEND_RUN_SLOW=1 (~45 s)")
    def test_base_e(self) -> None:
        # P(2) = e^^^e overflows, so no hex(2) check
        self._check("e", KN_REF_E, with_p2=False, ref_tol="1e-27")

    @unittest.skipUnless(RUN_SLOW, "Set FATOU_BACKEND_RUN_SLOW=1 (~25 s)")
    def test_base_10(self) -> None:
        """exp-081: centre and radius planned together. Around Re p the arc had
        to reach 0.85 Im p (K = 740, 4.3 min, 13.5 digits); now ~20 s, ~33 digits."""
        reg, err, h = _hex(self.fork, "10", ["hexreg", "hexerr", "hex(1/2)"])
        self.assertEqual(reg, 2)
        self.assertGreater(mp.re(err), 25)
        self.assertLess(abs(mp.re(h) - mp.mpf(KN_REF_10_HALF)), mp.mpf("1e-25"))


if __name__ == "__main__":
    unittest.main()
