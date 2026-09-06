from __future__ import annotations

import os
import sys
import unittest
from pathlib import Path

from mpmath import mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend import FatouGP, find_default_fatou_gp, find_default_gp_exe


@unittest.skipUnless(os.getenv("FATOU_BACKEND_RUN_SLOW") == "1", "Set FATOU_BACKEND_RUN_SLOW=1")
class FatouBackendSlowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.gp = FatouGP(
            gp_exe=find_default_gp_exe(),
            fatou_gp=find_default_fatou_gp(),
            dps=90,
            nlim=28,
            nskip=4,
            looplim=35,
        )
        mp.dps = 120

    def test_high_precision_complex_roundtrip_for_one_plus_i(self) -> None:
        y = mp.mpc("1.08", "0.12")
        residual = self.gp.roundtrip_residuals("1+I", [y])[0]
        self.assertLess(abs(residual), mp.mpf("1e-18"))

    def test_sub_eta_base_matches_independent_reference(self) -> None:
        """Regression guard for real bases below eta = e^(1/e).

        This used to compare the fork against the unmodified original. That
        oracle could not work: below eta there is no complex-conjugate fixed
        point for Kneser to use, and both engines returned the same value,
        correct to only ~15 digits, with no error raised
        (research/METHODS.md, section 2). Agreement between two engines that
        share a construction
        measures the construction, not the answer.

        The fork now uses regular iteration at the real attracting fixed
        point, so the oracle is the independent Koenigs/Schroeder value from
        research/tools/regular_subeta.py.
        """
        mp.dps = 60
        ref = mp.mpf("1.1362624867271280841853009186474260255893554650374189878532")
        fork = FatouGP(dps=38, fatou_gp="fork", state_cache=False)
        vf = fork.sexp("1.2", "0.5")
        self.assertLess(abs(vf - ref), mp.mpf("1e-35"))
        self.assertAlmostEqual(float(vf.real), 1.13626248673, places=9)
