from __future__ import annotations

import sys
import unittest
from pathlib import Path

from mpmath import mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend import FatouGP, find_default_fatou_gp, find_default_gp_exe
from fatou_backend.worker import FatouGPWorker, WorkerDied


def make_gp() -> FatouGP:
    return FatouGP(gp_exe=find_default_gp_exe(), fatou_gp=find_default_fatou_gp(),
                   dps=50, nlim=20, nskip=4, looplim=30)


class WorkerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        mp.dps = 80
        cls.gp = make_gp()
        cls.worker = FatouGPWorker(cls.gp.gp_exe, cls.gp._init_lines("e"))

    @classmethod
    def tearDownClass(cls) -> None:
        cls.worker.close()

    def test_eval_matches_oneshot(self) -> None:
        got = self.worker.eval(["sexp(0.5)"])[0]
        want = self.gp.sexp("e", mp.mpf("0.5"))
        self.assertLess(abs(got - want), mp.mpf("1e-30"))

    def test_second_eval_reuses_process(self) -> None:
        pid_before = self.worker._proc.pid
        values = self.worker.eval(["sexp(0.25)", "slog(2.0)"])
        self.assertEqual(len(values), 2)
        self.assertEqual(self.worker._proc.pid, pid_before)
        self.assertTrue(self.worker.alive)

    def test_gp_error_raises_with_output(self) -> None:
        with self.assertRaises(RuntimeError) as ctx:
            self.worker.eval(["undefined_function_xyz(1)"])
        self.assertIn("***", str(ctx.exception))
        # worker was killed on error; next use must fail fast
        self.assertFalse(self.worker.alive)
        with self.assertRaises(WorkerDied):
            self.worker.eval(["sexp(0.5)"])
        # respawn a shared worker for remaining tests in this class
        type(self).worker = FatouGPWorker(self.gp.gp_exe, self.gp._init_lines("e"))

    def test_timeout_kills_worker(self) -> None:
        w = FatouGPWorker(self.gp.gp_exe, self.gp._init_lines("e"), eval_timeout=3.0)
        with self.assertRaises(WorkerDied):
            w.eval(["while(1,)"])
        self.assertFalse(w.alive)
        w.close()

    def test_context_manager_closes(self) -> None:
        with FatouGPWorker(self.gp.gp_exe, self.gp._init_lines("e")) as w:
            w.eval(["sexp(0.1)"])
            self.assertTrue(w.alive)
        self.assertFalse(w.alive)


if __name__ == "__main__":
    unittest.main()
