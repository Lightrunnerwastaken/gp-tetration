"""Python and CLI access to pentation and hexation (exp-080).

The GP side is tested in test_pentation.py and test_hexation.py; this file
checks the plumbing: the backend methods return what the engine computes,
the per-process lazy initialisation survives a pool, the original engine is
refused loudly, and the CLI command works end to end.

Cost: ~20 s with a warm state cache.
"""
from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path

import mpmath as mp

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from fatou_backend import FatouGP, find_default_gp_exe

DPS = 30


class LevelsApiTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        mp.mp.dps = 80
        cls.fork = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=DPS)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.fork.close()

    def test_pent_matches_the_engine(self) -> None:
        api = self.fork.pent("2", mp.mpf("1.5"))
        direct = self.fork.eval_batch("2", ["pentinit(sexp(1),1)", "pent(3/2)"])[1]
        self.assertLess(abs(api - direct), mp.mpf(f"1e-{DPS}"))
        # Levenstein, Tetration Forum tid 1263 (2020): 2^^^1.5 ~ 2.7029
        self.assertAlmostEqual(float(mp.re(api)), 2.7029, places=4)

    def test_inverses_round_trip(self) -> None:
        x = mp.mpf("0.37")
        y = self.fork.pent("e", x)
        self.assertLess(abs(self.fork.ipent("e", y) - x), mp.mpf(f"1e-{DPS - 3}"))
        h = self.fork.hex("1.5", x)                   # regular hexation, b < eta5
        self.assertLess(abs(self.fork.ihex("1.5", h) - x), mp.mpf(f"1e-{DPS - 3}"))

    def test_hex_regular_identities(self) -> None:
        h0, h1, h2 = self.fork.hex_batch("1.5", [0, 1, 2])
        p15 = self.fork.pent("1.5", mp.mpf("1.5"))
        self.assertLess(abs(h0 - 1), mp.mpf(f"1e-{DPS - 3}"))
        self.assertLess(abs(h1 - mp.mpf("1.5")), mp.mpf(f"1e-{DPS - 3}"))
        self.assertLess(abs(h2 - p15), mp.mpf(f"1e-{DPS - 3}"))

    def test_digit_estimates(self) -> None:
        # the engine promises 0.95*dps; measured 29.1 (pentation, e) at dps 30
        self.assertGreater(self.fork.pentation_digits("e"), 0.9*DPS)
        self.assertGreater(self.fork.hexation_digits("1.5"), 0.9*DPS)

    def test_pool_initialises_each_worker(self) -> None:
        xs = [mp.mpf(k)/8 for k in range(-4, 8)]
        single = self.fork.pent_batch("2", xs)
        with FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=DPS, n_workers=2) as pool:
            pooled = pool.pent_batch("2", xs)
        for a, b in zip(single, pooled):
            self.assertLess(abs(a - b), mp.mpf(f"1e-{DPS}"))

    def test_original_engine_is_refused(self) -> None:
        orig = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="original", dps=DPS)
        try:
            with self.assertRaises(ValueError):
                orig.pent("2", mp.mpf("0.5"))
        finally:
            orig.close()

    def test_overflow_is_an_error_not_zero(self) -> None:
        # exp-083: at dps 16 the chain for 100^^^1.96875 ends in a real zero
        # without valid digits (0.E495855885102019914); Python read it as 0
        # and the calculator table carried that 0 behind its overflow gap
        with FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=16) as gp:
            with self.assertRaises(RuntimeError):
                gp.pent("100", mp.mpf("1.96875"))

    def test_cli(self) -> None:
        out = subprocess.run(
            [sys.executable, "-m", "fatou_backend.cli", "pent", "--base", "2",
             "--values", "1.5", "--fatou-gp", "fork", "--dps", str(DPS)],
            capture_output=True, text=True, cwd=REPO / "src", timeout=600)
        self.assertEqual(out.returncode, 0, out.stderr)
        self.assertTrue(out.stdout.strip().startswith("2.7028905714"), out.stdout)


if __name__ == "__main__":
    unittest.main()
