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


class CachedEqualsFreshTests(unittest.TestCase):
    """A restored engine must be indistinguishable from a freshly built one.

    Regression (2026-10-09): at dps 40 a base-e engine restored from the cache
    printed 40 digits where a fresh one printed 57 -- from bit-identical
    internal values. `default(realprecision, n)` also sets the printed digits,
    the fork's loop() ends with default(realprecision, precis), and \\uv does
    not see GP defaults. Runs at mpmath's default precision on purpose: the
    anchor check used to compare at mp.mp.dps, so at 15 the cache never hit
    and the restore path went untested.
    """
    DPS = 40
    EXPRS = ["sexp(1/2)", "slog(1/2)", "sexp(-17/10)", "slog(-185/100)",
             "sexp(1/3+I/5)"]
    PROBE = ("[default(realprecision), default(realbitprecision), "
             "default(seriesprecision), default(format)]")

    def setUp(self) -> None:
        self._tmp = tempfile.TemporaryDirectory()
        self._env = mock.patch.dict(os.environ, {"FATOU_CACHE_DIR": self._tmp.name})
        self._env.start()

    def tearDown(self) -> None:
        self._env.stop()
        self._tmp.cleanup()

    def _engine(self, cache: bool) -> FatouGP:
        return FatouGP(gp_exe=find_default_gp_exe(), fatou_gp="fork",
                       dps=self.DPS, state_cache=cache)

    def _run(self, base: str, cache: bool):
        with self._engine(cache) as gp:
            values = gp.eval_batch(base, self.EXPRS)
            worker = gp._workers[gp._base_expr(base)]
            probe = worker.eval_raw([f"print({self.PROBE})"], "__PROBE__", 60.0)
            restored = not any("sexpinit" in ln for ln in worker.debug_init_lines)
        return values, probe, restored

    def _assert_bit_identical(self, base: str) -> None:
        with mp.workdps(15):
            fresh, fresh_probe, _ = self._run(base, cache=False)
            self._run(base, cache=True)                       # fills the cache
            cached, cached_probe, restored = self._run(base, cache=True)
        self.assertTrue(restored, f"base {base}: second engine did not use the cache")
        self.assertEqual(cached_probe, fresh_probe)
        for expr, a, b in zip(self.EXPRS, fresh, cached):
            with self.subTest(base=base, value=expr):
                self.assertEqual(a, b, f"base {base}: {expr} differs, "
                                       f"rel {mp.nstr(abs(a - b) / abs(a), 3)}")

    def test_base_e(self) -> None:
        self._assert_bit_identical("e")

    def test_sub_eta_base(self) -> None:
        self._assert_bit_identical("1.2")


if __name__ == "__main__":
    unittest.main()
