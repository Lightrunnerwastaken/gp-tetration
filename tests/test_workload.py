from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from bench.workload import Case, all_cases, case_id


class WorkloadTests(unittest.TestCase):
    def test_default_cases_shape(self) -> None:
        cases = all_cases()
        tiers = {c.tier for c in cases}
        self.assertEqual(tiers, {"throughput", "precision", "anchor"})
        throughput = [c for c in cases if c.tier == "throughput"]
        self.assertEqual({c.base for c in throughput}, {"e", "2", "10"})
        self.assertTrue(all(c.dps == 80 for c in throughput))
        # 16 sexp- + 16 slog-Argumente pro Basis
        self.assertEqual(len(throughput), 3 * 32)
        precision = [c for c in cases if c.tier == "precision"]
        self.assertEqual({(c.base, c.dps) for c in precision}, {("e", 200), ("2", 200)})
        anchors = [c for c in cases if c.tier == "anchor"]
        self.assertEqual({c.base for c in anchors}, {"1+I", "0.8+0.4*I", "2+I"})

    def test_deep_adds_500_and_1000(self) -> None:
        deep = [c for c in all_cases(deep=True) if c.tier == "precision"]
        self.assertEqual({c.dps for c in deep}, {200, 500, 1000})

    def test_case_id_roundtrip_unique(self) -> None:
        cases = all_cases(deep=True)
        ids = [case_id(c) for c in cases]
        self.assertEqual(len(ids), len(set(ids)))


if __name__ == "__main__":
    unittest.main()
