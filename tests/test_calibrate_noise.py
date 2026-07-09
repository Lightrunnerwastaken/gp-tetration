from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "research"))

from calibrate_noise import noise_threshold_pct


class NoiseThresholdTests(unittest.TestCase):
    def test_floor_is_three_percent(self) -> None:
        self.assertEqual(noise_threshold_pct([1000.0, 1000.1, 999.9]), 3.0)

    def test_wide_spread_raises_threshold(self) -> None:
        # spread (1200-800)/1000 = 40% -> threshold 80%
        self.assertAlmostEqual(noise_threshold_pct([800.0, 1000.0, 1200.0]), 80.0)

    def test_requires_at_least_three_samples(self) -> None:
        with self.assertRaises(ValueError):
            noise_threshold_pct([1.0, 2.0])


if __name__ == "__main__":
    unittest.main()
