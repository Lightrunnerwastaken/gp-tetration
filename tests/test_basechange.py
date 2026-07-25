"""M5.1 tests: peel ladder Phi/mu against reference values published
in the mixed-base tetration paper series (see papers/)."""
import json
import os
import sys
import unittest
from pathlib import Path

import mpmath as mp

# This was the only test file without it, so it passed purely because the
# package happened to be installed editable.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend.gp_backend import FatouGP
from fatou_backend import basechange

MU_E2 = mp.mpf("-1.128403776628924009879217902036426267112")
MU_35 = mp.mpf("0.5788350021355712162374734720607948660945")


REPO = Path(__file__).resolve().parents[1]


def _proven(key: str, digits: int = 80) -> mp.mpf:
    """Load a proven reference without depending on the current directory.

    These two tests used a bare relative open("research/reference/..."), so they
    failed -- not skipped, failed -- whenever pytest ran from anywhere but the
    repo root. They are also the only two atlas tests that CAN fail, which made
    the failure mode particularly unhelpful.
    """
    payload = json.loads((REPO / "research" / "reference" / "values.json")
                         .read_text(encoding="utf-8"))
    return mp.mpf(payload["values"][key]["real"][:digits])


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
        ref = _proven("sexp|2|0.5|500")
        v = basechange.sexp_anchor(self.gp, 2, mp.mpf("0.5"))
        self.assertLess(abs(v - ref), mp.mpf("1e-20"))

    def test_sexp_anchor_matches_direct_engine(self):
        for y in ("0.0", "0.25", "0.9"):
            va = basechange.sexp_anchor(self.gp, 2, mp.mpf(y))
            vd = mp.mpf(self.gp.sexp(2, mp.mpf(y)).real)
            self.assertLess(abs(va - vd), mp.mpf("1e-20"))


    def test_slog_anchor_inverts_proven_value(self):
        ref = _proven("sexp|2|0.5|500")
        s = basechange.slog_anchor(self.gp, 2, ref)
        self.assertLess(abs(s - mp.mpf("0.5")), mp.mpf("1e-20"))

    def test_anchor_solver_inverts_itself(self):
        """Guards the Newton solve in slog_anchor -- NOT the atlas values.

        sexp_anchor adds Phi and slog_anchor solves it away with the SAME mode
        table, so any error in Phi cancels exactly. Demonstrated: perturbing mu
        for base 2 by 1e-9 moves the value-vs-reference check from 3.0e-24 to
        9.7e-10 while this residual stays at 3.8e-52, unchanged. The test that
        actually pins the atlas is test_sexp_anchor_matches_proven_reference.
        """
        for b in (2, 3):
            y = mp.mpf("0.3")
            rt = basechange.slog_anchor(self.gp, b,
                                        basechange.sexp_anchor(self.gp, b, y))
            self.assertLess(abs(rt - y), mp.mpf("1e-40"))
