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
  * by pinning what the builder must DELIVER (it must beat repeated
    multiplication, written out in the test rather than borrowed),
  * by asserting positively that every one of the seven twiddle sites calls the
    doubling builder, and
  * across two precisions, because the twiddle cache's precision key is
    load-bearing and single-precision tests cannot see it.

The earlier version of this file guarded the defect by grepping for the literal
`powers(`. That was defeated in review by a four-line local reimplementation
under another name: the shipped defect returned at its documented magnitude
with every test green. A name is not a property.

Runs in a few seconds; the one high-precision probe (n=1536 at precis 404) is
there because at dps 60 the effect being measured is ~1e-75 and the margin is
only 2.26x against a 2x threshold.
"""
from __future__ import annotations

import re
import subprocess
import sys
import unittest
from decimal import Decimal
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


def _decimals(out: str, tag: str) -> list[Decimal]:
    """Parse tagged GP values as Decimal.

    Required whenever the quantity can be smaller than ~1e-308: `float()` flushes
    those to 0.0, which turns a comparison of two error magnitudes into 0 < 0 and
    makes the test vacuous. At precis 404 the twiddle errors are ~1e-402, i.e.
    squarely in that region.
    """
    vals = []
    for line in out.splitlines():
        if line.startswith(tag + " "):
            vals.append(Decimal(line[len(tag) + 1:].strip()
                                .replace(" E", "E").replace(" ", "")))
    return vals


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

    def test_geopow_beats_repeated_multiplication(self) -> None:
        """Pin what the builder must DELIVER, not what it is called.

        Compared against repeated multiplication written out here, so the test
        states the property independently of PARI's `powers`: a doubling table
        carries ~log2(n) ulp where repeated multiplication carries ~n.

        Run at n=1536 / precis 404, where the measured advantage is ~19x. The
        earlier version used n=384 at dps 60, where the advantage is only 2.26x
        against a hardcoded 2.0 floor -- 12% of margin, one mantissa-granularity
        change away from flaking.
        """
        out = gp_run("\n".join([
            "mulpow(mu, n) = {",
            "  local(v, j); v = vector(n+1); v[1] = 1;",
            "  for (j = 1, n, v[j+1] = v[j] * mu); return(v);",
            "};",
            "n = 1536;",
            # Reference roots at a precision well above the working one --
            # `precision(Pi, 900)` would only PAD a 404-digit Pi, not recompute it.
            "default(realprecision, 520);",
            "refv = vector(n, j, exp(-2*Pi*I*(j-1)/n));",
            "default(realprecision, 404);",
            "mu = exp(-2*Pi*I/n);",
            "a = mulpow(mu, n-1); b = geopow(mu, n-1);",
            "ep = 0; eg = 0;",
            # One line on purpose: gp reading a script from stdin does NOT
            # continue an open parenthesis across a newline unless it is inside
            # braces. Split, this loop never runs and the probe reports 0/0 --
            # which is why assertGreater(ep, 0) below is load-bearing, not decor.
            "for (j = 1, n, ep = max(ep, abs(a[j]-refv[j])); eg = max(eg, abs(b[j]-refv[j])));",
            'print("EP ", ep); print("EG ", eg);',
        ]), dps=404)
        ep, eg = _decimals(out, "EP"), _decimals(out, "EG")
        self.assertTrue(ep and eg, f"probe produced no output:\n{out}")
        self.assertGreater(ep[0], 0, "the probe measured no error at all -- it "
                                     "did not run; see _decimals on why float() "
                                     "cannot be used here")
        self.assertLess(
            eg[0], ep[0] / 5,
            f"the twiddle builder is no better than repeated multiplication "
            f"(repeated {ep[0]:.3e}, builder {eg[0]:.3e}) -- the doubling "
            f"construction has been lost")

    def test_geoseq_is_c_times_powers(self) -> None:
        out = gp_run(
            "c = exp(0.3*I); mu = exp(2*Pi*I/61); n = 50;"
            "v = geoseq(c, mu, n);"
            'print("LEN ", if(#v == n, 1, 0));'
            'print("DEV ", vecmax(vector(n, j, abs(v[j] - c*mu^j))));'
        )
        self.assertEqual(_floats(out, "LEN"), [1.0], f"bad length:\n{out}")
        self.assertLess(_floats(out, "DEV")[0], 1e-50)


def _strip_gp_comments(src: str) -> str:
    """Blank comment bodies but KEEP their newlines, so line numbers survive.

    Deleting `/* ... */` outright shifts every subsequent line number (measured:
    a 97-line drift by the middle of the fork), which makes any reported
    offender location wrong. PARI's other comment form is `\\\\ ...` to
    end-of-line.
    """
    src = re.sub(r"/\*.*?\*/", lambda m: "\n" * m.group().count("\n"),
                 src, flags=re.DOTALL)
    return re.sub(r"\\\\.*", "", src)


class TwiddleSiteTests(unittest.TestCase):
    """Guard the twiddle construction at each SITE, positively.

    The first version of this guard grepped for the literal `powers(` and was
    trivially defeatable: `powers()` is not a language feature, it is repeated
    multiplication, so a four-line local reimplementation under any other name
    reintroduces the defect with the guard green. That was demonstrated against
    this file -- the shipped defect came back at its documented magnitude
    (34x worse at n=4608) with all tests passing.

    So the guard is now positive: each of the seven sites must be calling the
    doubling builder. Swapping in any other construction removes the expected
    call and fails here. Paired with `test_geopow_beats_repeated_multiplication`
    below, which pins what the builder must actually DELIVER, the two together
    cover both "is it wired up" and "is it correct".
    """

    # (description, exact source fragment that must be present)
    SITES = [
        ("thtaylor FFT roots", "G = fft(geopow(om^(-1), samples-1), t_est)"),
        ("thtaylor output map", "pw = geopow(conj(c0)/om, terms);"),
        ("staylor extraction, complex", "G = fft(geopow(om^(-1), samples-1), exin)"),
        ("staylor extraction, real", "G = fft(geopow(mu^(-1), 2*samples-1), concat(exin"),
        ("mixdft inner FFT roots", "w = geopow(exp(-2*Pi*I/m), m-1);"),
        ("mixdft outer twiddles", "geopow(exp(-2*Pi*I/nn), nn-1)"),
        ("bluedft FFT roots", "w = geopow(exp(2*Pi*I/M), M-1);"),
    ]

    def test_every_twiddle_site_uses_the_doubling_builder(self) -> None:
        src = _strip_gp_comments(FORK.read_text(encoding="utf-8", errors="replace"))
        missing = [name for name, frag in self.SITES if frag not in src]
        self.assertEqual(
            missing, [],
            "these twiddle sites no longer build their tables with geopow():\n  "
            + "\n  ".join(missing)
            + "\nRepeated multiplication accumulates ~n ulp; the doubling build "
              "carries ~log2(n). Renaming the helper does not make it correct.",
        )

    def test_no_live_powers_call_in_the_fork(self) -> None:
        """Cheap supplement to the positive check above, with honest line numbers."""
        stripped = _strip_gp_comments(FORK.read_text(encoding="utf-8", errors="replace"))
        offenders = [f"line {i}: {ln.strip()}"
                     for i, ln in enumerate(stripped.splitlines(), 1)
                     if "powers(" in ln]
        self.assertEqual(offenders, [],
                         "fatou_fork.gp calls powers() outside a comment:\n"
                         + "\n".join(offenders))

    def test_comment_stripping_preserves_line_numbers(self) -> None:
        src = FORK.read_text(encoding="utf-8", errors="replace")
        self.assertEqual(len(_strip_gp_comments(src).splitlines()),
                         len(src.splitlines()),
                         "comment stripping shifted line numbers, so every "
                         "offender location this file reports would be wrong")


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

    def test_twiddle_cache_respects_precision(self) -> None:
        """The cache's precision key is load-bearing; nothing else covered it.

        mixdft/bluedft cache their twiddle tables on (nn, realprecision), and
        the engine ladders precision *within* a run -- the extraction transforms
        at a reduced `exdig` while the theta side runs at full precision. If the
        precision half of the key is dropped, a table built cheap gets reused
        expensive and the result silently loses everything past the low
        precision. Measured cost of removing that clause: ~232 digits.

        Every other test here runs at a single precision, so this was invisible.
        """
        out = gp_run("\n".join([
            self.NAIVE,
            # prime the cache at LOW precision, then use the same nn high
            "default(realprecision, 60);",
            "ta = tvec(96); xa = mixdft(ta);",
            "default(realprecision, 300);",
            "tb = tvec(96); xb = mixdft(tb); xc = naive(tb);",
            'print("HI ", vecmax(vector(96, j, abs(xb[j]-xc[j]))));',
            # and a genuine same-nn cache HIT at one precision
            "td = tvec(96); xd = mixdft(td); xe = naive(td);",
            'print("HIT ", vecmax(vector(96, j, abs(xd[j]-xe[j]))));',
        ]), dps=60)
        hi, hit = _floats(out, "HI"), _floats(out, "HIT")
        self.assertTrue(hi and hit, f"probe produced no output:\n{out}")
        self.assertLess(hi[0], 1e-250,
                        f"a twiddle table built at low precision was reused at "
                        f"high precision (error {hi[0]:.3e}) -- the precision "
                        f"clause of the cache key has been lost")
        self.assertLess(hit[0], 1e-250,
                        f"repeat call on a cached grid is wrong: {hit[0]:.3e}")

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
