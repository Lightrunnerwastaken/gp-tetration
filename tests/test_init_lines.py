from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend import FatouGP, find_default_fatou_gp, find_default_gp_exe


class InitLinesTests(unittest.TestCase):
    def setUp(self) -> None:
        self.gp = FatouGP(
            gp_exe=find_default_gp_exe(),
            fatou_gp=find_default_fatou_gp(),
            dps=50,
            nlim=20,
            nskip=4,
            looplim=30,
        )

    def test_init_lines_order_and_content(self) -> None:
        lines = self.gp._init_lines("e")
        self.assertEqual(lines[0], "default(realprecision, 50);")
        self.assertTrue(lines[1].startswith('read("'))
        self.assertTrue(lines[1].endswith('");'))
        self.assertIn("fatou.gp", lines[1])
        self.assertNotIn("\\", lines[1])  # as_posix, keine Backslashes
        self.assertEqual(lines[2], "quietmode=1;")
        self.assertEqual(lines[3], "sexpinit(exp(1),20,4,30);")
        self.assertEqual(len(lines), 4)

    def test_init_lines_complex_base(self) -> None:
        lines = self.gp._init_lines("1+I")
        self.assertEqual(lines[3], "sexpinit(1+I,20,4,30);")


if __name__ == "__main__":
    unittest.main()
