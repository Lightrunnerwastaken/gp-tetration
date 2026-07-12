"""M5.1 tests: peel ladder Phi/mu against values published by the
basechange-modes research log (read-only reference values)."""
import unittest

import mpmath as mp

from fatou_backend.gp_backend import FatouGP
from fatou_backend import basechange

MU_E2 = mp.mpf("-1.128403776628924009879217902036426267112")
MU_35 = mp.mpf("0.5788350021355712162374734720607948660945")


def _gp() -> FatouGP:
    return FatouGP(dps=60)


class BasechangeMuTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 70
        cls.gp = _gp()

    def test_mu_e2_matches_research_log(self):
        m = basechange.mu(self.gp, self.gp, "exp(1)", 2, n_grid=16)
        self.assertLess(abs(m - MU_E2), mp.mpf("1e-14"))

    def test_mu_35_matches_research_log(self):
        m = basechange.mu(self.gp, self.gp, 3, 5, n_grid=16)
        self.assertLess(abs(m - MU_35), mp.mpf("1e-14"))

    def test_mu_antisymmetry_e2(self):
        m_fwd = basechange.mu(self.gp, self.gp, "exp(1)", 2, n_grid=16)
        m_rev = basechange.mu(self.gp, self.gp, 2, "exp(1)", n_grid=16)
        self.assertLess(abs(m_fwd + m_rev), mp.mpf("1e-14"))


if __name__ == "__main__":
    unittest.main()
