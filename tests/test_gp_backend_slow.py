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

    def test_sub_eta_base_matches_original_engine(self) -> None:
        """Regression: exp-032 once broke real bases below eta = e^(1/e).

        The gate never covers that regime, so guard it here via engine
        diversity: fork and unmodified original must agree at b=1.2
        (attracting-fixed-point tetration, sexp(0.5) ~ 1.13626).
        """
        mp.dps = 60
        fork = FatouGP(dps=38, fatou_gp="fork", state_cache=False)
        orig = FatouGP(dps=38, fatou_gp="original", state_cache=False)
        vf = fork.sexp("1.2", "0.5")
        vo = orig.sexp("1.2", "0.5")
        self.assertLess(abs(vf - vo), mp.mpf("1e-25"))
        self.assertAlmostEqual(float(vf.real), 1.13626248673, places=9)
