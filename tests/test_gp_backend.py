from __future__ import annotations

import os
import sys
import unittest
from pathlib import Path
from unittest import mock

from mpmath import mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend import FatouGP, find_default_fatou_gp, find_default_gp_exe

FIXTURE_FATOU_GP = Path(__file__).resolve().parent / "fixtures" / "override_fatou.gp"


class FatouBackendTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.gp = FatouGP(
            gp_exe=find_default_gp_exe(),
            fatou_gp=find_default_fatou_gp(),
            dps=50,
            nlim=20,
            nskip=4,
            looplim=30,
        )
        mp.dps = 80

    @classmethod
    def tearDownClass(cls) -> None:
        cls.gp.close()

    def test_sexp_e_half_matches_known_value(self) -> None:
        value = self.gp.sexp("e", mp.mpf("0.5"))
        self.assertLess(abs(value - mp.mpf("1.64635423375119458097")), mp.mpf("1e-18"))

    def test_slog_sexp_roundtrip_for_real_value(self) -> None:
        x = mp.mpf("0.35")
        y = self.gp.sexp("2", x)
        recovered = self.gp.slog("2", y)
        # Measured residual at these settings (dps 50, nlim 20, looplim 30):
        # 2.69e-43. The old 1e-18 left 25 digits of slack -- an engine that
        # lost 25 digits passed. Threshold is measured-worst x 100.
        self.assertLess(abs(recovered - x), mp.mpf("1e-40"))

    def test_complex_values_parse(self) -> None:
        value = self.gp.sexp("e", mp.mpc("0.5", "0.2"))
        self.assertTrue(mp.isfinite(mp.re(value)))
        self.assertTrue(mp.isfinite(mp.im(value)))

    def test_vendored_fatou_gp_is_default(self) -> None:
        with mock.patch.dict(os.environ, {}, clear=False):
            path = find_default_fatou_gp()
        self.assertEqual(path.name, "fatou.gp")
        self.assertIn("vendor", path.parts)

    def test_env_override_for_fatou_gp_takes_priority(self) -> None:
        with mock.patch.dict(os.environ, {"FATOU_GP_FILE": str(FIXTURE_FATOU_GP)}, clear=False):
            self.assertEqual(find_default_fatou_gp(), FIXTURE_FATOU_GP)

    def test_explicit_paths_take_priority(self) -> None:
        gp = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp=FIXTURE_FATOU_GP, dps=40)
        self.assertEqual(gp.fatou_gp, FIXTURE_FATOU_GP)

    def test_missing_gp_exe_raises_clear_error(self) -> None:
        with self.assertRaises(FileNotFoundError):
            FatouGP(gp_exe=Path(r"C:\definitely_missing\gp.exe"), fatou_gp=find_default_fatou_gp())

    def test_missing_fatou_gp_raises_clear_error(self) -> None:
        with self.assertRaises(FileNotFoundError):
            FatouGP(gp_exe=find_default_gp_exe(), fatou_gp=Path(r"C:\definitely_missing\fatou.gp"))

    def test_session_helper_matches_direct_calls(self) -> None:
        session = self.gp.session("1+I")
        batch = session.sexp_batch([mp.mpf("0.25"), mp.mpc("0.4", "0.2")])
        direct = [
            self.gp.sexp("1+I", mp.mpf("0.25")),
            self.gp.sexp("1+I", mp.mpc("0.4", "0.2")),
        ]
        for left, right in zip(batch, direct):
            self.assertLess(abs(left - right), mp.mpf("1e-20"))

    def test_complex_eval_batch_parser_handles_genuine_complex_output(self) -> None:
        value = self.gp.eval_batch("1+I", ["sexp(0.5)"])[0]
        self.assertTrue(mp.isfinite(mp.re(value)))
        self.assertTrue(mp.isfinite(mp.im(value)))
        self.assertGreater(abs(mp.im(value)), mp.mpf("1e-12"))

    def test_complex_batch_matches_single_for_multiple_bases(self) -> None:
        bases = ["1+I", "0.8+0.4*I", "2+I"]
        xs = [mp.mpf("0.5"), mp.mpc("0.4", "0.2")]
        for base in bases:
            batch = self.gp.sexp_batch(base, xs)
            single = [self.gp.sexp(base, x) for x in xs]
            for left, right in zip(batch, single):
                self.assertLess(abs(left - right), mp.mpf("1e-18"))
                self.assertTrue(mp.isfinite(mp.re(left)))
                self.assertTrue(mp.isfinite(mp.im(left)))

    def test_complex_roundtrip_residuals_are_small(self) -> None:
        cases = {
            "1+I": [mp.mpc("1.05", "0.1"), mp.mpc("1.2", "-0.2")],
            "0.8+0.4*I": [mp.mpc("0.95", "0.15")],
            "2+I": [mp.mpc("1.1", "0.05")],
        }
        for base, ys in cases.items():
            residuals = self.gp.roundtrip_residuals(base, ys)
            for residual in residuals:
                # Measured worst residual across these four cases: 8.26e-52
                # (1+I); the others reach 1e-56/1e-57. The old 1e-12 left 40-45
                # digits of slack, on the only numeric assertion the default
                # suite makes about complex bases. Measured-worst x 100.
                self.assertLess(abs(residual), mp.mpf("1e-48"))


if __name__ == "__main__":
    unittest.main()
