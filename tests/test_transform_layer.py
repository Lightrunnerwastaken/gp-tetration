"""Direct tests for the fork's transform layer (geopow/geoseq/mixdft/bluedft).

Why these exist. The 2026-07-25 review found a ~2.7-digit accuracy loss in the
DFT twiddle factors that had shipped as a "verified" keep. Nothing caught it:

  * the end-to-end gate passed -- a few ulp in the twiddles do not move the
    final sexp value far enough to trip a threshold calibrated on digits;
  * the keep's own verification compared mixdft against bluedft, and *both*
    built their twiddles with `powers()`, so the two agreed on a common-mode
    error and the agreement was read as proof.

So the transform layer had no test that could see a wrong answer at the level
where it was wrong. These tests check it directly and independently:

  * against a naive O(n^2) DFT (a different algorithm, not a shared helper),
  * against high-precision truth for the power tables,
  * and, for the twiddles, they assert that the doubling build is measurably
    BETTER than `powers()` -- which is the regression guard: if someone puts
    `powers()` back, this fails.

Kept small (dps 60, n <= 384) so the whole file runs in a few seconds.
"""
from __future__ import annotations

import re
import subprocess
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fatou_backend import find_default_gp_exe

REPO = Path(__file__).resolve().parents[1]
FORK = REPO / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"


def gp_run(body: str, dps: int = 60) -> str:
    """Run a GP snippet against the fork and return stdout."""
    script = "\n".join([
        "default(parisizemax, 2147483648);",
        f"default(realprecision, {dps});",
        f'read("{FORK.as_posix()}");',
        f"default(realprecision, {dps});",
        "quietmode=1;",
        body,
        "quit;",
    ]) + "\n"
    proc = subprocess.run(
        [str(find_default_gp_exe()), "-q", "-f"],
        input=script, capture_output=True, text=True,
        encoding="utf-8", errors="replace", cwd=str(REPO),
    )
    return proc.stdout


def _floats(out: str, tag: str) -> list[float]:
    vals = []
    for line in out.splitlines():
        if line.startswith(tag + " "):
            vals.append(float(line[len(tag) + 1:].strip()
                              .replace(" E", "e").replace(" ", "")))
    return vals


class PowerTableTests(unittest.TestCase):
    """geopow must reproduce powers() semantics and beat its accuracy."""

    def test_geopow_matches_powers_semantics(self) -> None:
        # powers(mu, n) is [1, mu, ..., mu^n]; geopow must be the same length
        # and the same values to working precision.
        out = gp_run(
            "mu = exp(2*Pi*I/97);"
            "a = powers(mu, 40); b = geopow(mu, 40);"
            'print("LEN ", if(#a == #b, 1, 0));'
            'print("DEV ", vecmax(vector(#a, j, abs(a[j]-b[j]))));'
        )
        self.assertEqual(_floats(out, "LEN"), [1.0], f"length mismatch:\n{out}")
        dev = _floats(out, "DEV")
        self.assertTrue(dev, f"no DEV line:\n{out}")
        self.assertLess(dev[0], 1e-50, "geopow disagrees with powers()")

    def test_geopow_is_more_accurate_than_powers(self) -> None:
        """The regression guard for the 2026-07-25 twiddle defect.

        powers() accumulates ~n ulp, a doubling table ~log2(n). If someone
        swaps geopow back to powers at a twiddle site, this is what fails.
        """
        out = gp_run(
            "n = 384; mu = exp(-2*Pi*I/n);"
            "a = powers(mu, n-1); b = geopow(mu, n-1);"
            "ep = 0; eg = 0;"
            "for (j = 1, n,"
            "  ref = exp(precision(-2*Pi*I, 200)*(j-1)/n);"
            "  ep = max(ep, abs(a[j]-ref)); eg = max(eg, abs(b[j]-ref));"
            ");"
            'print("EP ", ep); print("EG ", eg);'
        )
        ep, eg = _floats(out, "EP"), _floats(out, "EG")
        self.assertTrue(ep and eg, f"probe produced no output:\n{out}")
        self.assertGreater(ep[0], 0.0)
        self.assertLess(eg[0], ep[0] / 2.0,
                        "geopow is not measurably better than powers() -- "
                        "a twiddle site has regressed to powers()")

    def test_geoseq_is_c_times_powers(self) -> None:
        out = gp_run(
            "c = exp(0.3*I); mu = exp(2*Pi*I/61); n = 50;"
            "v = geoseq(c, mu, n);"
            'print("LEN ", if(#v == n, 1, 0));'
            'print("DEV ", vecmax(vector(n, j, abs(v[j] - c*mu^j))));'
        )
        self.assertEqual(_floats(out, "LEN"), [1.0], f"bad length:\n{out}")
        self.assertLess(_floats(out, "DEV")[0], 1e-50)


class TwiddleSiteTests(unittest.TestCase):
    """The actual regression guard: no live `powers()` left in the fork.

    Worth being explicit about why this is a SOURCE check and not a numeric
    one. The DFT tests below run at dps 60, where a `powers()` twiddle table of
    length 288 is off by ~288 ulp ~ 3e-65 -- three orders of magnitude below
    their 1e-45 tolerance. They cannot see this defect, and raising the
    precision far enough to see it would make them slow. What they do catch is
    a gross algorithmic error; what catches the twiddle defect is this.

    Seven sites used `powers()` on 2026-07-25 and every one of them fed either
    an FFT root table or an output coefficient mapping. `geopow` is the
    doubling replacement; if any site regresses, this fails.
    """

    def test_no_live_powers_call_in_the_fork(self) -> None:
        src = FORK.read_text(encoding="utf-8", errors="replace")
        stripped = re.sub(r"/\*.*?\*/", "", src, flags=re.DOTALL)
        offenders = [
            f"line {i}: {ln.strip()}"
            for i, ln in enumerate(stripped.splitlines(), 1)
            if "powers(" in ln
        ]
        self.assertEqual(
            offenders, [],
            "fatou_fork.gp calls powers() outside a comment. powers() "
            "accumulates ~n ulp in the twiddle tables; use geopow(), which "
            "builds the same vector by doubling (~log2(n) ulp). Offenders:\n"
            + "\n".join(offenders),
        )

    def test_geopow_is_defined(self) -> None:
        src = FORK.read_text(encoding="utf-8", errors="replace")
        self.assertIn("geopow(mu, n) = {", src)


class ExactDFTTests(unittest.TestCase):
    """mixdft/bluedft against a naive DFT -- a genuinely different algorithm.

    Convention (from the fork): X[k] = sum_j t[j] * omega^((j-1)(k-1)),
    omega = exp(-2*Pi*I/N).
    """

    # Joined with newlines and terminated with ';' -- concatenated without a
    # separator, "}tvec(...)" is a syntax error and gp answers with its banner.
    NAIVE = "\n".join([
        "naive(t) = {",
        "  local(nn, om);",
        "  nn = #t; om = exp(-2*Pi*I/nn);",
        # reduce the exponent mod nn: keeps the powering cheap and avoids
        # comparing against a reference that is itself sloppy.
        "  return(vector(nn, kk, sum(j = 1, nn,",
        "    t[j] * om^(((j-1)*(kk-1)) % nn))));",
        "};",
        "tvec(nn) = vector(nn, j, (j % 7) - 3 + I*((j % 5) - 2));",
        "",
    ])

    def _worst(self, fn: str, sizes: list[int]) -> dict[int, float]:
        # Hoist the transforms out of the vector body -- written inline they
        # would be re-evaluated once per element, i.e. n DFTs instead of one.
        body = self.NAIVE + "\n".join(
            f"ta = tvec({n}); xa = {fn}(ta); xb = naive(ta);\n"
            f'print("D{n} ", vecmax(vector({n}, j, abs(xa[j] - xb[j]))));'
            for n in sizes
        )
        out = gp_run(body)
        got = {}
        for n in sizes:
            vals = _floats(out, f"D{n}")
            self.assertTrue(vals, f"no output for n={n}:\n{out}")
            got[n] = vals[0]
        return got

    def test_mixdft_matches_naive_dft(self) -> None:
        # 96 = 3*32, 160 = 5*32, 288 = 9*32 exercise the radix-r path;
        # 128 is r == 1 (pure power of two); 176 = 11*16 has r > 9 and must
        # fall back to bluedft. All must be exact to working precision.
        for n, err in self._worst("mixdft", [96, 128, 160, 176, 288]).items():
            self.assertLess(err, 1e-45, f"mixdft wrong at n={n}: {err}")

    def test_bluedft_matches_naive_dft(self) -> None:
        for n, err in self._worst("bluedft", [30, 96, 100]).items():
            self.assertLess(err, 1e-45, f"bluedft wrong at n={n}: {err}")

    def test_mixdft_and_bluedft_agree(self) -> None:
        # The original (insufficient) check, kept because it is still a useful
        # cross-check -- but it can only ever be a supplement to the tests
        # above, since both routines share the twiddle construction.
        out = gp_run(
            self.NAIVE +
            "ta = tvec(288); xa = mixdft(ta); xb = bluedft(ta);\n"
            'print("DEV ", vecmax(vector(288, j, abs(xa[j] - xb[j]))));'
        )
        self.assertLess(_floats(out, "DEV")[0], 1e-45)


if __name__ == "__main__":
    unittest.main()
