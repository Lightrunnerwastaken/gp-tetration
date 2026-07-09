from __future__ import annotations

import sys
import unittest
from pathlib import Path

from mpmath import mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend import FatouGP, find_default_fatou_gp, find_default_gp_exe


def make_gp(**kwargs) -> FatouGP:
    return FatouGP(gp_exe=find_default_gp_exe(), fatou_gp=find_default_fatou_gp(),
                   dps=50, nlim=20, nskip=4, looplim=30, **kwargs)


class PersistentIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        mp.dps = 80

    def test_persistent_matches_oneshot_real_base(self) -> None:
        with make_gp() as pers:
            oneshot = make_gp(persistent=False)
            xs = [mp.mpf("0.25"), mp.mpf("0.5")]
            a = pers.sexp_batch("e", xs)
            b = oneshot.sexp_batch("e", xs)
            for left, right in zip(a, b):
                self.assertLess(abs(left - right), mp.mpf("1e-30"))

    def test_persistent_matches_oneshot_complex_base(self) -> None:
        with make_gp() as pers:
            oneshot = make_gp(persistent=False)
            a = pers.sexp("1+I", mp.mpc("0.4", "0.2"))
            b = oneshot.sexp("1+I", mp.mpc("0.4", "0.2"))
            self.assertLess(abs(a - b), mp.mpf("1e-30"))

    def test_worker_is_reused_across_batches(self) -> None:
        with make_gp() as gp:
            gp.sexp("e", mp.mpf("0.25"))
            worker1 = gp._workers[gp._base_expr("e")]
            gp.slog("e", mp.mpf("2.0"))
            worker2 = gp._workers[gp._base_expr("e")]
            self.assertIs(worker1, worker2)

    def test_respawn_after_worker_death(self) -> None:
        with make_gp() as gp:
            gp.sexp("e", mp.mpf("0.25"))
            gp._workers[gp._base_expr("e")].close()  # simulate crash
            value = gp.sexp("e", mp.mpf("0.5"))
            self.assertLess(abs(value - mp.mpf("1.64635423375119458097")),
                            mp.mpf("1e-18"))

    def test_close_kills_all_workers(self) -> None:
        gp = make_gp()
        gp.sexp("e", mp.mpf("0.25"))
        workers = list(gp._workers.values())
        gp.close()
        self.assertTrue(all(not w.alive for w in workers))
        self.assertEqual(gp._workers, {})


if __name__ == "__main__":
    unittest.main()
