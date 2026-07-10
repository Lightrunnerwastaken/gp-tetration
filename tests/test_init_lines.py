from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend import FatouGP, find_default_fatou_gp, find_default_gp_exe
from fatou_backend.gp_backend import GP_EXE_CANDIDATES


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
        self.assertEqual(lines[0], "default(parisizemax, 2147483648);")
        self.assertEqual(lines[1], "default(realprecision, 50);")
        self.assertTrue(lines[2].startswith('read("'))
        self.assertTrue(lines[2].endswith('");'))
        self.assertIn("fatou.gp", lines[2])
        self.assertNotIn("\\", lines[2])  # as_posix, keine Backslashes
        self.assertEqual(lines[3], "quietmode=1;")
        self.assertEqual(lines[4], "sexpinit(exp(1),20,4,30);")
        self.assertEqual(len(lines), 5)

    def test_init_lines_complex_base(self) -> None:
        lines = self.gp._init_lines("1+I")
        self.assertEqual(lines[4], "sexpinit(1+I,20,4,30);")


class GpExeCandidateTests(unittest.TestCase):
    def test_pari64_is_preferred(self) -> None:
        self.assertIn("Pari64", str(GP_EXE_CANDIDATES[0]))
        self.assertIn("Pari64", str(GP_EXE_CANDIDATES[1]))
        self.assertIn("Pari32", str(GP_EXE_CANDIDATES[2]))

    def test_default_gp_exe_resolves_to_pari64(self) -> None:
        self.assertIn("Pari64", str(find_default_gp_exe()))


if __name__ == "__main__":
    unittest.main()
