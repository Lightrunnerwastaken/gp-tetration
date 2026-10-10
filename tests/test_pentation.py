"""Pentation in the fork (exp-076): precision, base coverage, independent checks.

The original fatou.gp builds the Poincare function of slog at the negative
fixed point of sexp with a FIXED 21-term series, so its pentation carries
~30 digits whatever precision is asked for (dps 60: pent(0)-1 = 6.5e-31), and
it cannot run below eta at all -- its error estimate compares with `z1<1`
and slogtaylor is guarded there. The fork scales the series with the
working precision and goes through slog(), which routes to the regular
iteration for 1 < b < eta.

What is checked, and against what:

  * the defining identities -- P(0)=1, P(1)=b, P(-1)=0, P(-2)=-1,
    P(2)=sexp(b) and the inverse roundtrip -- on real bases on both sides of
    eta, a large base and a complex base. They test the series against each
    other and against the engine, not the engine's own digits;
  * the error vector: the same values at dps 40 and dps 60. Different working
    precision, different slog samples, different series degree -- agreement
    bounds the error of the shallower run (research/METHODS.md S1);
  * an independent method: J. D. Fox (Tetration Forum tid 103, 2007) computed
    L, lambda and e^^^n with a 1200-term matrix slog, ~22 correct digits;
  * engine diversity: the unmodified original's pentation at its ~30 digits.

Cost: ~20 s cold, a few seconds with the state cache.
"""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend import FatouGP, find_default_gp_exe

DPS = 40
DEEP = 60

# pentinit(B, 1): B is only used to decide whether the series may be taken
# real, so the engine's own sexp(1) = b serves for every base expression.
INIT = "pentinit(sexp(1),1)"

# Fox 2007, as quoted in the tetration-graph experiment exp.pentation. The
# trailing digits are his; ~22 of them are correct.
FOX = {
    "L": "-1.850354529027181418483437788",
    "lambda": "6.460671295681839390208370083",
    "P(2)": "2075.968335058065833574141757",
    "P(-3)": "-1.636358354286028979629049436",
    "P(-4)": "-1.813170483098635639971748853",
    "P(-10)": "-1.850354005584711078364293582",
}


def _pent(engine: FatouGP, base: str, exprs: list[str]) -> list[mp.mpc]:
    """pentinit in the same worker, then the expressions; drops errz."""
    return engine.eval_batch(base, [INIT] + exprs)[1:]


class PentationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        mp.mp.dps = 120
        cls.fork = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=DPS)
        cls.deep = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=DEEP)
        cls.orig = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="original", dps=DPS)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.fork.close()
        cls.deep.close()
        cls.orig.close()

    # -- identities ---------------------------------------------------------
    # Measured at dps 40 (working precision 57 digits): the worst of the
    # residuals below is 7.5e-48 (base e, P(0)-1, fresh engine), most sit at
    # 1e-55..1e-57. They test internal consistency; a cache-restored engine
    # gives bit-identical values (see FLOOR below).
    # The threshold is the engine's own promise (~dps digits), not the
    # residuals: the pentation cannot be better than the sexp it iterates.
    TOL = mp.mpf(f"1e-{DPS}")

    def _assert_identities(self, base: str, with_p2: bool = True, tol=None) -> None:
        tol = self.TOL if tol is None else tol
        exprs = ["pent(0)-1", "pent(1)-sexp(1)", "pent(-1)", "pent(-2)+1",
                 "ipent(pent(1/3))-1/3"]
        if with_p2:
            exprs.append("pent(2)/sexp(sexp(1))-1")
        for expr, val in zip(exprs, _pent(self.fork, base, exprs)):
            with self.subTest(base=base, check=expr):
                self.assertLess(abs(val), tol,
                                f"base {base}: {expr} = {mp.nstr(val, 5)}")

    def test_identities_base_e(self) -> None:
        self._assert_identities("e")

    def test_identities_base_2(self) -> None:
        self._assert_identities("2")

    def test_identities_sub_eta(self) -> None:
        """1 < b < eta: the original cannot run here at all."""
        self._assert_identities("1.2")

    def test_identities_large_base(self) -> None:
        """Base 100: lambda = 1100, L only 2e-4 above the singularity at -2.

        Regression guard. The Newton solve leaves the high Poincare
        coefficients with 19 of 57 digits -- legitimately, their
        contribution |c|*pentr^n is below 1e-65 -- but PARI carries the
        precision through Horner RELATIVELY, and pent(0)-1 came out 3.6e-21.
        pentclean pads such terms only when their error contribution is
        provably below the target. P(2) = 100^^100 is not representable.

        Measured: 3.0e-57 (P(0)-1), worst 1.6e-51 (roundtrip). Without the
        padding 2.6e-40, without padding and truncation 3.6e-21 -- hence a
        threshold well below the general one.
        """
        self._assert_identities("100", with_p2=False, tol=mp.mpf("1e-48"))

    def test_identities_complex_base(self) -> None:
        self._assert_identities("2+I/2")

    def test_error_estimate_is_reported(self) -> None:
        # pentinit returns its own digit estimate (the original returned 0
        # whenever sexpinit was skipped). Measured 47.1 for base e at dps 40.
        errz = mp.re(self.fork.eval_batch("e", [INIT])[0])
        self.assertGreater(errz, DPS, f"pentinit estimates only {errz} digits")

    # -- error vector: dps 40 against dps 60 ----------------------------------
    # The threshold is the engine's documented floor (0.95*dps digits). Until
    # 2026-10-09 it also had to absorb the state cache: a restored engine
    # printed only dps digits (40-41 against a fresh dps-100 run, base e) where
    # a fresh one printed its working precision (54-56), because the cache did
    # not carry GP's `format`. Fixed; cached and fresh engines are now
    # bit-identical (tests/test_state_cache.py, CachedEqualsFreshTests), so the
    # digits here no longer depend on whether the cache is warm. The pentation
    # inherits whatever the engine delivers.
    FLOOR = mp.mpf(10) ** (-0.95 * DPS)

    def test_error_vector_against_deeper_run(self) -> None:
        exprs = ["pent(1/2)", "pent(-1/2)", "pent(3/2)", "tfixed", "lambdat"]
        for base in ("e", "1.2"):
            shallow = _pent(self.fork, base, exprs)
            deep = _pent(self.deep, base, exprs)
            for expr, a, b in zip(exprs, shallow, deep):
                with self.subTest(base=base, value=expr):
                    rel = abs(a - b) / abs(b)
                    self.assertLess(rel, self.FLOOR,
                                    f"base {base}: {expr} moves by {mp.nstr(rel, 5)} "
                                    f"between dps {DPS} and {DEEP}")

    # -- independent method -------------------------------------------------
    def test_matches_fox_2007(self) -> None:
        names = ["L", "lambda", "P(2)", "P(-3)", "P(-4)", "P(-10)"]
        exprs = ["tfixed", "lambdat", "pent(2)", "pent(-3)", "pent(-4)", "pent(-10)"]
        for name, val in zip(names, _pent(self.fork, "e", exprs)):
            with self.subTest(value=name):
                ref = mp.mpf(FOX[name])
                rel = abs(mp.re(val) - ref) / abs(ref)
                self.assertLess(rel, mp.mpf("1e-20"),
                                f"{name}: {mp.nstr(rel, 5)} away from Fox (2007)")

    def test_matches_levenstein_base_2(self) -> None:
        # Tetration Forum tid 1263 (2020): 2^^^1.5 ~ 2.7029, L ~ -1.7439.
        p15, L = _pent(self.fork, "2", ["pent(3/2)", "tfixed"])
        self.assertAlmostEqual(float(mp.re(p15)), 2.7029, places=4)
        self.assertAlmostEqual(float(mp.re(L)), -1.7439, places=4)

    # -- engine diversity ---------------------------------------------------
    def test_fork_matches_original_pentation(self) -> None:
        """The original's fixed 21-term series carries ~30 digits at dps 40."""
        exprs = ["pent(1/2)", "pent(-1/2)", "tfixed"]
        for base in ("e", "2"):
            new = _pent(self.fork, base, exprs)
            old = _pent(self.orig, base, exprs)
            for expr, a, b in zip(exprs, new, old):
                with self.subTest(base=base, value=expr):
                    self.assertLess(abs(a - b), mp.mpf("1e-25"),
                                    f"base {base}: fork and original disagree on {expr}")


if __name__ == "__main__":
    unittest.main()
