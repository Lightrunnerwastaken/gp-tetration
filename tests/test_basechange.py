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


class PhiModeTableTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 70
        cls.gp = _gp()

    def test_reconstruction_matches_ladder(self):
        import tempfile, os
        path = os.path.join(tempfile.gettempdir(), "phi_modes_test.json")
        if os.path.exists(path):
            os.remove(path)
        mu_v, modes = basechange.phi_modes_cached(self.gp, 2, path=path)
        for th in ("0.037", "0.41", "0.777"):
            direct = basechange.phi(self.gp, self.gp, "exp(1)", 2, mp.mpf(th))
            recon = basechange.phi_from_modes(mu_v, modes, mp.mpf(th))
            self.assertLess(abs(direct - recon), mp.mpf("1e-20"))
        # Cache-Hit liefert identische Werte
        mu_c, modes_c = basechange.phi_modes_cached(self.gp, 2, path=path)
        self.assertLess(abs(mu_c - mu_v), mp.mpf("1e-40"))
        os.remove(path)


class SexpAnchorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        mp.mp.dps = 70
        cls.gp = _gp()

    def test_sexp_anchor_matches_proven_reference(self):
        import json
        ref = mp.mpf(json.load(open("research/reference/values.json"))
                     ["values"]["sexp|2|0.5|500"]["real"][:80])
        v = basechange.sexp_anchor(self.gp, 2, mp.mpf("0.5"))
        self.assertLess(abs(v - ref), mp.mpf("1e-20"))

    def test_sexp_anchor_matches_direct_engine(self):
        for y in ("0.0", "0.25", "0.9"):
            va = basechange.sexp_anchor(self.gp, 2, mp.mpf(y))
            vd = mp.mpf(self.gp.sexp(2, mp.mpf(y)).real)
            self.assertLess(abs(va - vd), mp.mpf("1e-20"))
