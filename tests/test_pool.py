from __future__ import annotations

import sys
import unittest
from pathlib import Path

from mpmath import mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend import FatouGP, find_default_fatou_gp, find_default_gp_exe
from fatou_backend.pool import FatouGPPool


def make_gp(**kwargs) -> FatouGP:
    return FatouGP(gp_exe=find_default_gp_exe(), fatou_gp=find_default_fatou_gp(),
                   dps=50, nlim=20, nskip=4, looplim=30, **kwargs)


class PoolTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        mp.dps = 80

    def test_pool_preserves_order_and_values(self) -> None:
        with make_gp() as single:
            xs = [mp.mpf(j) / 8 for j in range(8)]
            want = single.sexp_batch("e", xs)
        with make_gp(n_workers=2) as pooled:
            got = pooled.sexp_batch("e", xs)
        self.assertEqual(len(got), len(want))
        for left, right in zip(got, want):
            self.assertLess(abs(left - right), mp.mpf("1e-30"))

    def test_pool_class_round_robin(self) -> None:
        gp = make_gp()
        pool = FatouGPPool(lambda: gp._spawn_worker("e"), n_workers=2)
        try:
            values = pool.eval([f"sexp({j}/8)" for j in range(5)])
            self.assertEqual(len(values), 5)
            direct = gp._spawn_worker("e")
            try:
                expected = direct.eval([f"sexp({j}/8)" for j in range(5)])
            finally:
                direct.close()
            for left, right in zip(values, expected):
                self.assertLess(abs(left - right), mp.mpf("1e-30"))
        finally:
            pool.close()
            gp.close()

    def test_small_batch_avoids_pool(self) -> None:
        with make_gp(n_workers=4) as gp:
            gp.sexp("e", mp.mpf("0.5"))   # 1 expression < 2*4 -> single worker
            self.assertEqual(gp._pools, {})


if __name__ == "__main__":
    unittest.main()
