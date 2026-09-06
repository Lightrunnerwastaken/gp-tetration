"""End-to-end coverage of the OPTIMIZED fork -- which the suite had none of.

Every other test helper builds its engine with `find_default_fatou_gp()`, and
that resolves to `vendor/fatou.gp`, the pristine original. The three files that
name the fork at all are `test_transform_layer.py` (unit-level: it never calls
sexpinit/sexp/slog), `test_gate.py` and `test_gp_backend_slow.py` -- and the
latter two are skipped unless FATOU_BACKEND_RUN_SLOW=1, which nothing in the
repo sets. So the artifact that carries 30 optimizations plus today's
correctness fix had zero end-to-end coverage in the suite that gates the work.

What this file adds:

  * engine diversity -- fork against the unmodified original at identical
    settings. Two independent code paths agreeing to ~50 digits is the
    strongest practical evidence that the optimizations did not bend accuracy,
    and it is the same argument METHODS makes for the reference ladder.
  * the sub-eta regime (1 < b < e^(1/e)), which the exp-032/033/034 cache
    extensions once silently broke and which is otherwise covered only by a
    skipped test.
  * the published calibration claim (>= 0.95*dps true digits), measured against
    the proven references rather than asserted in prose.

Cost: ~30 s the first time (cold engine init per base/engine), ~0.3 s
afterwards -- the state cache makes repeats essentially free.
"""
from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend import FatouGP, find_default_gp_exe

REPO = Path(__file__).resolve().parents[1]
REFERENCE = REPO / "research" / "reference" / "values.json"
DPS = 50


def _reference(key: str, digits: int = 70) -> mp.mpf:
    """A proven reference value, truncated well inside its proven depth."""
    payload = json.loads(REFERENCE.read_text(encoding="utf-8"))
    return mp.mpf(payload["values"][key]["real"][:digits])


class ForkEngineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        mp.mp.dps = 120
        cls.fork = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=DPS)
        cls.orig = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="original", dps=DPS)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.fork.close()
        cls.orig.close()

    # -- engine diversity ---------------------------------------------------
    # Measured agreement at dps 50: 4.2e-50 (e), 7.1e-51 (2). The threshold
    # sits ~4 orders above the worst of those: tight enough that a real
    # regression (which collapses agreement by tens of digits) fails, loose
    # enough not to flake on a different PARI build.
    # b=1.2 used to be listed here at 3.6e-51. That number was real but
    # meaningless: below eta both engines run the same construction and agree
    # on its error. See test_fork_matches_independent_reference_sub_eta.
    DIVERSITY_TOL = mp.mpf("1e-45")

    def _assert_engines_agree(self, base: str, x: str) -> None:
        a = self.fork.sexp(base, mp.mpf(x))
        b = self.orig.sexp(base, mp.mpf(x))
        delta = abs(a - b)
        self.assertLess(
            delta, self.DIVERSITY_TOL,
            f"fork and original disagree for base {base} at x={x}: "
            f"{mp.nstr(delta, 6)} -- an optimization has bent accuracy")

    def test_fork_matches_original_base_e(self) -> None:
        self._assert_engines_agree("e", "0.5")

    def test_fork_matches_original_base_2(self) -> None:
        self._assert_engines_agree("2", "0.5")

    # sexp_1.2(0.5) by regular iteration at the real attracting fixed point
    # (Koenigs/Schroeder, research/tools/regular_subeta.py; self-tests S(1)=b
    # and S(2)=b^b hold to 148 digits -- S(0)=1 is tautological here and is
    # deliberately not used as evidence). Owes nothing to fatou.gp.
    # Kept as a STRING on purpose: a class-level mp.mpf(...) is evaluated at
    # import time, when mp.mp.dps is still the default 15, and the constant
    # would silently arrive rounded to 15 digits.
    SUB_ETA_REF_STR = (
        "1.136262486727128084185300918647426025589355465037418987853243230891553944148904")

    def test_fork_matches_independent_reference_sub_eta(self) -> None:
        """1 < b < e^(1/e): the attracting-fixed-point regime.

        This test used to assert fork == original, and that oracle was wrong.
        Kneser needs a complex-conjugate fixed-point pair, which does not
        exist below eta, and BOTH engines returned the same value that is
        only ~15 digits correct -- silently, with an imaginary part of
        7.26e-14 on a provably real quantity (research/METHODS.md S2). Two
        implementations of one construction agree on their shared error,
        which is exactly what engine diversity is meant to exclude and
        cannot, when the two engines share the construction.

        The fork now routes these bases through regular iteration, so the
        oracle here is an independent method, and the original is kept only
        for the shallow cross-check it can actually support.
        """
        ref = mp.mpf(self.SUB_ETA_REF_STR)
        got = self.fork.sexp("1.2", mp.mpf("0.5"))
        self.assertLess(
            abs(got - ref), mp.mpf(f"1e-{DPS - 3}"),
            "the fork lost the sub-eta regime")
        self.assertLess(
            abs(mp.im(mp.mpc(got))), mp.mpf(f"1e-{DPS - 3}"),
            "sexp is real for 1 < b < eta; an imaginary part is the defect's tell")

        old = self.orig.sexp("1.2", mp.mpf("0.5"))
        self.assertLess(
            abs(old - ref), mp.mpf("1e-12"),
            "the original should still be right in its leading digits")
        self.assertGreater(
            abs(old - ref), mp.mpf(f"1e-{DPS - 3}"),
            "the original's sub-eta defect is gone -- did upstream fix it? "
            "If so this test and METHODS section 2 both need revisiting")

    def test_fork_matches_original_complex_base(self) -> None:
        self._assert_engines_agree("1+I", "0.5")

    # sexp_1.4494(0.5). Generated on the fork after the exp-074 fix at dps 80
    # (90.6 contour digits); the leading ~60 are independently confirmed by the
    # unmodified original, which reaches 60.67 digits there and agrees to all
    # of them. String, not mp.mpf -- see SUB_ETA_REF_STR above for why.
    NEAR_ETA_REF_STR = (
        "1.2592243909108569064268681775693833372881536039326069230")

    def test_near_eta_base_survives(self) -> None:
        """b just above eta = 1.4446678610: the near-boundary band.

        Two keeps have silently broken a base regime the gate does not cover.
        exp-032 broke the sub-eta side (fixed 2026-07-15); exp-048 broke this
        one, found 2026-07-30 by sweeping all 34 fork versions. Its degree
        truncation anchored the tolerance to the largest *coefficient*, which
        is only the largest *contribution* when the coefficients decay -- i.e.
        when circr > 1. Base e has circr = 1.3372, this base has ~0.133, so
        the coefficients grow, mx sat at the top index, and the threshold came
        out ~85 decimal digits too high: 1.837 digits instead of 147.4 at
        dps 150, with no error raised.

        The gate cannot catch this: all six of its bases sit comfortably far
        from eta. So it is guarded here, and deliberately with a value rather
        than a digit count -- a wrong value is the symptom, a low digit count
        only its cause.
        """
        ref = mp.mpf(self.NEAR_ETA_REF_STR)
        got = self.fork.sexp("1.4494", mp.mpf("0.5"))
        self.assertLess(
            abs(got - ref), mp.mpf(f"1e-{DPS - 5}"),
            "the near-eta band regressed -- check icbuild's truncation "
            "threshold (exp-074) before anything else")

    # -- the published calibration claim ------------------------------------
    def test_fork_delivers_the_claimed_true_digits(self) -> None:
        """README/CHANGELOG promise >= 0.95*dps true digits. Measure it.

        Measured at dps 50: 53.7 digits for base e, 53.1 for base 2 -- the
        engine overshoots its nominal dps, which is why the claim is a floor.
        Asserted against the proven references, not against the other engine.
        """
        floor = 0.95 * DPS
        for base, key in (("e", "sexp|e|0.5|500"), ("2", "sexp|2|0.5|500")):
            with self.subTest(base=base):
                got = mp.re(self.fork.sexp(base, mp.mpf("0.5")))
                err = abs(got - _reference(key))
                digits = float(-mp.log10(err)) if err else float("inf")
                self.assertGreaterEqual(
                    digits, floor,
                    f"base {base}: {digits:.1f} true digits at dps {DPS}, "
                    f"below the documented floor of {floor:.1f}")

    # -- the fork must round-trip -------------------------------------------
    def test_fork_roundtrip_real_and_complex(self) -> None:
        # Measured residuals at these settings: 2.4e-42 (base 2, real),
        # 3.5e-42 / 3.0e-43 / 9.1e-43 (complex bases). Threshold two orders up.
        for base in ("2", "1+I", "2+I", "0.8+0.4*I"):
            with self.subTest(base=base):
                x = mp.mpf("0.35") if base == "2" else mp.mpf("0.5")
                back = self.fork.slog(base, self.fork.sexp(base, x))
                self.assertLess(abs(back - x), mp.mpf("1e-40"),
                                f"fork roundtrip broken for base {base}")


if __name__ == "__main__":
    unittest.main()
