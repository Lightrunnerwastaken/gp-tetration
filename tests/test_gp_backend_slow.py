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
