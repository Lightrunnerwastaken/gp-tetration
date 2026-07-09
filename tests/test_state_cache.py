from __future__ import annotations

import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from mpmath import mp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend import FatouGP, find_default_fatou_gp, find_default_gp_exe
from fatou_backend import state_cache


def make_gp(**kwargs) -> FatouGP:
    return FatouGP(gp_exe=find_default_gp_exe(), fatou_gp=find_default_fatou_gp(),
                   dps=50, nlim=20, nskip=4, looplim=30, **kwargs)


class StateCacheTests(unittest.TestCase):
    def setUp(self) -> None:
        mp.dps = 80
        self._tmp = tempfile.TemporaryDirectory()
        self._env = mock.patch.dict(os.environ, {"FATOU_CACHE_DIR": self._tmp.name})
        self._env.start()

    def tearDown(self) -> None:
        self._env.stop()
        self._tmp.cleanup()

    def test_cache_key_depends_on_dps(self) -> None:
        gp = make_gp()
        k1 = state_cache.cache_key("exp(1)", 50, 20, 4, 30, gp.fatou_gp, gp.gp_exe)
        k2 = state_cache.cache_key("exp(1)", 60, 20, 4, 30, gp.fatou_gp, gp.gp_exe)
        self.assertNotEqual(k1, k2)
        self.assertEqual(len(k1), 64)

    def test_cache_written_on_first_use_and_used_on_second(self) -> None:
        with make_gp() as gp1:
            v1 = gp1.sexp("e", mp.mpf("0.5"))
        bins = list(Path(self._tmp.name).glob("*.gpbin"))
        self.assertEqual(len(bins), 1)
        self.assertTrue(bins[0].with_suffix(".json").exists())
        with make_gp() as gp2:
            v2 = gp2.sexp("e", mp.mpf("0.5"))
            # cached worker must not have run sexpinit
            worker = gp2._workers[gp2._base_expr("e")]
            self.assertFalse(any("sexpinit" in ln for ln in worker.debug_init_lines))
        self.assertLess(abs(v1 - v2), mp.mpf("1e-30"))

    def test_corrupt_cache_falls_back_to_full_init(self) -> None:
        with make_gp() as gp1:
            v1 = gp1.sexp("e", mp.mpf("0.5"))
        bin_path = next(Path(self._tmp.name).glob("*.gpbin"))
        bin_path.write_bytes(b"garbage")
        with make_gp() as gp2:
            v2 = gp2.sexp("e", mp.mpf("0.5"))
        self.assertLess(abs(v1 - v2), mp.mpf("1e-30"))

    def test_state_cache_disabled(self) -> None:
        with make_gp(state_cache=False) as gp:
            gp.sexp("e", mp.mpf("0.5"))
        self.assertEqual(list(Path(self._tmp.name).glob("*")), [])


if __name__ == "__main__":
    unittest.main()
