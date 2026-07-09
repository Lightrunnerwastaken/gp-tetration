from __future__ import annotations

import sys
import unittest
from pathlib import Path

import mpmath as mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from bench.metrics import correct_digits


class MetricsTests(unittest.TestCase):
    def setUp(self) -> None:
        mp.mp.dps = 120

    def test_exact_match_returns_cap(self) -> None:
        v = mp.mpc("1.5", "0.25")
        self.assertEqual(correct_digits(v, v, cap=80), 80.0)

    def test_known_relative_error(self) -> None:
        ref = mp.mpc("2.0", "0")
        val = ref + mp.mpf("2e-10")  # rel err 1e-10
        d = correct_digits(val, ref, cap=80)
        self.assertAlmostEqual(d, 10.0, places=3)

    def test_small_reference_uses_absolute_scale(self) -> None:
        ref = mp.mpc("1e-30", "0")   # |ref| < 1 -> scale = 1
        val = ref + mp.mpf("1e-12")
        d = correct_digits(val, ref, cap=80)
        self.assertAlmostEqual(d, 12.0, places=3)

    def test_cap_is_upper_bound(self) -> None:
        ref = mp.mpc("2.0", "0")
        val = ref + mp.mpf("1e-300")
        self.assertEqual(correct_digits(val, ref, cap=50), 50.0)


if __name__ == "__main__":
    unittest.main()
