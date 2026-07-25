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
    # Measured agreement at dps 50: 4.2e-50 (e), 7.1e-51 (2), 3.6e-51 (1.2).
    # The threshold sits ~4 orders above the worst of those: tight enough that
    # a real regression (which collapses agreement by tens of digits) fails,
    # loose enough not to flake on a different PARI build.
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

    def test_fork_matches_original_sub_eta_base(self) -> None:
        """1 < b < e^(1/e): the attracting-fixed-point regime.

        This is the regime the exp-032/033/034 cache extensions broke without
        the gate noticing (the gate does not test it), which is why the fork
        routes these bases through the pre-extension path via the `subeta`
        flag. Nothing in the default suite exercised that flag.
        """
        self._assert_engines_agree("1.2", "0.5")

    def test_fork_matches_original_complex_base(self) -> None:
        self._assert_engines_agree("1+I", "0.5")

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
