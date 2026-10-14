"""Hyperoperations of any integer rank (exp-088).

H_0 = x+1, H_1 = b+x, H_2 = b*x, H_3 = b^x, H_4 = sexp, H_5 = pentation,
H_6 = hexation; from rank 7 on the fork builds H_r from H_{r-1}: odd ranks
regular at the left real fixed point of H_{r-1}, even ranks regular at its
attracting fixed point (b < eta_{r-1}); above eta_{r-1} even ranks are
Kneser at the complex fixed point of H_{r-1} (exp-089).

Reference values from the tetration-graph session
(forschung/eta_hyperoperationen, 12.10.2026), computed independently on a
snapshot of the fork with their own hyperoperation code:
  H_7(2, 0.5) = 1.48726615560199, H_7(2, 1.5) = 2.65161799652261
  H_8(2, 0.5) = 1.49025612402668, H_9(2, 0.5) = 1.49281110382266
  H_10(1.83, 0.5) ~ 1.4567116, H_10(1.83, 1.5) ~ 2.0842796 (backward Koenigs
  iteration, ~8 digits)

Cost: ~40 s (base 2 builds its Kneser hexation first); octation and
nonation at base 2 and decation at base 1.83 only with
FATOU_BACKEND_RUN_SLOW=1 (several minutes).
"""
from __future__ import annotations

import os
import subprocess
import sys
import unittest
from pathlib import Path

import mpmath as mp

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from fatou_backend import FatouGP, find_default_gp_exe

DPS = 20
RUN_SLOW = os.getenv("FATOU_BACKEND_RUN_SLOW") == "1"


class HyperTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        mp.mp.dps = 60
        cls.fork = FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=DPS)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.fork.close()

    def test_elementary_ranks(self) -> None:
        h = lambda r, x: self.fork.hyper("2", r, mp.mpf(x))
        self.assertLess(abs(h(0, "0.5") - mp.mpf("1.5")), mp.mpf("1e-30"))
        self.assertLess(abs(h(1, "0.5") - mp.mpf("2.5")), mp.mpf("1e-30"))
        self.assertLess(abs(h(2, "1.5") - 3), mp.mpf("1e-30"))
        self.assertLess(abs(h(3, "0.5") - mp.sqrt(2)), mp.mpf("1e-30"))
        self.assertLess(abs(self.fork.ihyper("2", 3, 8) - 3), mp.mpf("1e-30"))

    def test_ranks_4_to_6_are_the_existing_levels(self) -> None:
        x = mp.mpf("0.37")
        self.assertLess(abs(self.fork.hyper("2", 4, x) - self.fork.sexp("2", x)), mp.mpf("1e-30"))
        self.assertLess(abs(self.fork.hyper("2", 5, x) - self.fork.pent("2", x)), mp.mpf("1e-30"))
        self.assertLess(abs(self.fork.hyper("2", 6, x) - self.fork.hex("2", x)), mp.mpf("1e-28"))

    def test_heptation_base_2(self) -> None:
        # measured at dps 20: 29.6 digits, identities ~1e-31, init 0.3 s
        d = self.fork.hyper_digits("2", 7)
        self.assertGreater(d, 25)
        tol = mp.mpf(10) ** (-int(d) + 2)
        h0, h1, h2, h05, h15 = self.fork.hyper_batch("2", 7, [0, 1, 2, mp.mpf("0.5"), mp.mpf("1.5")])
        self.assertLess(abs(h0 - 1), tol)
        self.assertLess(abs(h1 - 2), tol)
        self.assertLess(abs(h2 - 4), tol)            # H_r(2, 2) = 4 for every r >= 2
        self.assertLess(abs(h05 - mp.mpf("1.48726615560199")), mp.mpf("1e-13"))
        self.assertLess(abs(h15 - mp.mpf("2.65161799652261")), mp.mpf("1e-13"))
        self.assertLess(abs(self.fork.ihyper("2", 7, h05) - mp.mpf("0.5")), tol)

    def test_regular_ranks_base_1_5(self) -> None:
        # b = 1.5 lies below every eta_r: all ranks regular. Each one is built
        # from the one below, so the recursion H_r(x+1) = H_{r-1}(H_r(x)) is
        # an independent check. Measured: ~32 digits up to rank 8.
        for r in (7, 8):
            with self.subTest(rank=r):
                x = mp.mpf("0.3")
                a, b = self.fork.hyper_batch("1.5", r, [x, x + 1])
                low = self.fork.hyper("1.5", r - 1, a)
                self.assertLess(abs(b - low), mp.mpf("1e-25"))
                self.assertGreater(self.fork.hyper_digits("1.5", r), 25)

    def test_sub_eta_base(self) -> None:
        # b = 6/5 < e^(1/e): tetration itself regular (diversity probe, see
        # STATE.md). Measured at dps 20: rank 7 32.7, rank 8 31.4 digits.
        for r in (7, 8):
            with self.subTest(rank=r):
                h0, h1 = self.fork.hyper_batch("6/5", r, [0, 1])
                self.assertLess(abs(h0 - 1), mp.mpf("1e-25"))
                self.assertLess(abs(h1 - mp.mpf(6)/5), mp.mpf("1e-25"))
                self.assertGreater(self.fork.hyper_digits("6/5", r), 25)

    def test_grid_matches_direct(self) -> None:
        # exp-091: hypergridval (calculator tables) takes one step of H_{r-1}
        # from the point one period earlier instead of a full chain per point
        step = mp.mpf(1)/8
        exprs = [f"hypergridval(7, -3., {mp.nstr(step, 5)}, {i})" for i in range(1, 41)]
        grid = self.fork.eval_batch("1.5", exprs)
        for i in (5, 12, 20, 33, 40):
            with self.subTest(i=i):
                direct = self.fork.hyper("1.5", 7, -3 + (i - 1)*step)
                self.assertLess(abs(grid[i - 1] - direct), mp.mpf("1e-25"))

    def test_complex_base_stops_at_hexation(self) -> None:
        # complex bases have only Sheldon's ~1-digit hexation; the fork
        # refuses to build rank 7 on it instead of returning noise
        with self.assertRaises(RuntimeError):
            self.fork.hyper("2+I/2", 7, mp.mpf("0.5"))

    @unittest.skipUnless(RUN_SLOW, "Set FATOU_BACKEND_RUN_SLOW=1 (~3-4 min)")
    def test_octation_and_nonation_base_2(self) -> None:
        # exp-089: b = 2 > eta_7 = 1.8398, so octation is Kneser at the complex
        # fixed point p = -4.2948 + 1.2835i of heptation -- the same
        # construction as the hexation on the pentation, now with a state per
        # rank. Nonation is regular on top of it. Measured at dps 20: octation
        # 29.5 digits. The graph session's values carry ~14-16 digits.
        refs = {8: ("1.49025612402668", "2.63586067475121"),
                9: ("1.49281110382266", "2.62453210644146")}
        for r, (r05, r15) in refs.items():
            with self.subTest(rank=r):
                d = self.fork.hyper_digits("2", r)
                self.assertGreater(d, 20)
                h0, h1, h2, h05, h15 = self.fork.hyper_batch(
                    "2", r, [0, 1, 2, mp.mpf("0.5"), mp.mpf("1.5")])
                tol = mp.mpf(10) ** (-int(d) + 2)
                self.assertLess(abs(h0 - 1), tol)
                self.assertLess(abs(h1 - 2), tol)
                self.assertLess(abs(h2 - 4), tol)
                self.assertLess(abs(h05 - mp.mpf(r05)), mp.mpf("1e-13"))
                self.assertLess(abs(h15 - mp.mpf(r15)), mp.mpf("1e-13"))

    def test_rank_must_be_a_nonnegative_integer(self) -> None:
        for bad in (-1, 2.5, True):
            with self.subTest(rank=bad):
                with self.assertRaises(ValueError):
                    self.fork.hyper("2", bad, mp.mpf("0.5"))

    @unittest.skipUnless(RUN_SLOW, "Set FATOU_BACKEND_RUN_SLOW=1 (~1.5 min)")
    def test_decation_base_2(self) -> None:
        # exp-089: Kneser on nonation. Its second arc fails (points near the
        # complex fixed point of octation cannot be inverted there) and a
        # shifted sibling arc certifies it. Reference: the same build at
        # dps 20 (28.2 digits); at dps 16 the build carries ~9.7 digits.
        with FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=16) as gp:
            self.assertGreater(gp.hyper_digits("2", 10), 8)
            h05, h15, h2 = gp.hyper_batch("2", 10, [mp.mpf("0.5"), mp.mpf("1.5"), 2])
        self.assertLess(abs(h05 - mp.mpf("1.4961951114859126262")), mp.mpf("1e-8"))
        self.assertLess(abs(h15 - mp.mpf("2.6186419110477304859")), mp.mpf("1e-8"))
        self.assertLess(abs(h2 - 4), mp.mpf("1e-8"))

    @unittest.skipUnless(RUN_SLOW, "Set FATOU_BACKEND_RUN_SLOW=1 (~1-2 min)")
    def test_rank_12_base_2(self) -> None:
        # exp-090: the Kneser arc of rank 12 lies around Re ~ -7; normalising
        # from 1 needs ~8 inverse steps, the old limit of 6 (tuned for the
        # hexation) refused it. At dps 16 the build carries ~8 digits.
        with FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork", dps=16) as gp:
            self.assertGreater(gp.hyper_digits("2", 12), 6)
            h0, h1, h2 = gp.hyper_batch("2", 12, [0, 1, 2])
        for got, want in ((h0, 1), (h1, 2), (h2, 4)):
            self.assertLess(abs(got - want), mp.mpf("1e-6"))

    @unittest.skipUnless(RUN_SLOW, "Set FATOU_BACKEND_RUN_SLOW=1 (~1.5 min)")
    def test_decation_base_1_83(self) -> None:
        # b = 1.83 < eta_7: octation regular, decation regular at the
        # attracting fixed point of nonation -- the same construction as the
        # graph's backward Koenigs iteration, which carried ~8 digits.
        h05, h15 = self.fork.hyper_batch("1.83", 10, [mp.mpf("0.5"), mp.mpf("1.5")])
        self.assertLess(abs(h05 - mp.mpf("1.4567116")), mp.mpf("1e-7"))
        self.assertLess(abs(h15 - mp.mpf("2.0842796")), mp.mpf("1e-7"))
        self.assertGreater(self.fork.hyper_digits("1.83", 10), 25)

    def test_cli(self) -> None:
        out = subprocess.run(
            [sys.executable, "-m", "fatou_backend.cli", "hyper", "--rank", "3", "--base", "2",
             "--values", "0.5", "--fatou-gp", "fork", "--dps", str(DPS)],
            capture_output=True, text=True, cwd=REPO / "src", timeout=600)
        self.assertEqual(out.returncode, 0, out.stderr)
        self.assertTrue(out.stdout.strip().startswith("1.414213562373"), out.stdout)


if __name__ == "__main__":
    unittest.main()
