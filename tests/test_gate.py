from __future__ import annotations

import os
import shutil
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "research"))

from gate import run_gate

REPO = Path(__file__).resolve().parents[1]
ORIGINAL = REPO / "src" / "fatou_backend" / "vendor" / "fatou.gp"
FORK = REPO / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"

RUN_SLOW = os.getenv("FATOU_BACKEND_RUN_SLOW") == "1"


class GateTests(unittest.TestCase):
    def test_fork_starts_byte_identical(self) -> None:
        self.assertEqual(ORIGINAL.read_bytes(), FORK.read_bytes())

    @unittest.skipUnless(RUN_SLOW, "Set FATOU_BACKEND_RUN_SLOW=1")
    def test_gate_passes_on_unmodified_fork(self) -> None:
        ok, report = run_gate(FORK)
        self.assertTrue(ok, "\n".join(report))

    @unittest.skipUnless(RUN_SLOW, "Set FATOU_BACKEND_RUN_SLOW=1")
    def test_gate_fails_on_sabotaged_fork(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sabotaged = Path(tmp) / "fatou_sabotaged.gp"
            shutil.copy(FORK, sabotaged)
            with sabotaged.open("a", encoding="utf-8") as handle:
                handle.write("\nslog(z) = { abel(z*lnb+k-1)+rslog + 1e-20; }\n")
            ok, report = run_gate(sabotaged)
            self.assertFalse(ok)
            self.assertTrue(any("FAIL" in line for line in report))


if __name__ == "__main__":
    unittest.main()
