from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Sequence

import mpmath as mp

from . import state_cache
from .pool import FatouGPPool
from .worker import FatouGPWorker, WorkerDied


GPValue = mp.mpf | mp.mpc | float | complex | int | str


def _first_existing(candidates: Iterable[Path]) -> Path | None:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


GP_EXE_CANDIDATES = (
    Path(r"C:\Program Files\Pari64-2-17-3\gp.exe"),
    # the NSIS installer is a 32-bit stub, so its default target is the x86 dir
    Path(r"C:\Program Files (x86)\Pari64-2-17-3\gp.exe"),
    Path(r"C:\Program Files (x86)\Pari32-2-17-3\gp.exe"),
)


def find_default_gp_exe() -> Path:
    # A set-but-wrong env var is a typo, not a request to fall back: silently
    # ignoring it means you run a different binary than you asked for.
    env = os.getenv("FATOU_GP_EXE")
    if env:
        path = Path(env)
        if not path.exists():
            raise FileNotFoundError(
                f"FATOU_GP_EXE is set to {env!r}, which does not exist. "
                f"Fix it or unset it; refusing to silently use another binary."
            )
        return path
    path = _first_existing(GP_EXE_CANDIDATES)
    if path is not None:
        return path
    # PARI/GP installs to PATH as plain `gp` on Linux/macOS (and via MSYS on
    # Windows); the hardcoded candidates above only cover the Windows installer.
    found = shutil.which("gp") or shutil.which("gp.exe")
    if found:
        return Path(found)
    raise FileNotFoundError(
        "Could not locate a PARI/GP executable. Tried $FATOU_GP_EXE, the "
        "standard Windows install paths, and `gp` on PATH. Install PARI/GP "
        "(apt install pari-gp / brew install pari / windows installer) or set "
        "FATOU_GP_EXE to the executable."
    )


def find_default_fatou_gp() -> Path:
    # Critical: the fallback below is the PRISTINE fatou.gp. Silently ignoring a
    # bad FATOU_GP_FILE means a typo makes you benchmark the unmodified original
    # while believing you are measuring the fork -- the one mistake this repo
    # can least afford.
    env = os.getenv("FATOU_GP_FILE")
    if env:
        path = Path(env)
        if not path.exists():
            raise FileNotFoundError(
                f"FATOU_GP_FILE is set to {env!r}, which does not exist. "
                f"Fix it or unset it; refusing to silently fall back to the "
                f"unmodified vendored engine."
            )
        return path
    here = Path(__file__).resolve()
    candidates = (
        here.parent / "vendor" / "fatou.gp",
        here.parents[2] / "vendor" / "fatou.gp",
    )
    path = _first_existing(candidates)
    if path is None:
        raise FileNotFoundError(
            "Could not locate fatou.gp. Set FATOU_GP_FILE or use the vendored file in the repo."
        )
    return path


def _parse_gp_scalar(text: str) -> mp.mpf:
    """Parse a GP-printed real at a precision that cannot truncate it.

    mpmath parses at the GLOBAL `mp.mp.dps`, which defaults to 15. Parsing a
    300-digit GP result under that default silently yields a double: the caller
    pays for the full high-precision engine run and gets ~16 true digits back,
    with everything past digit 17 being the binary expansion of that double --
    plausible-looking, and wrong. The string itself states how much precision it
    needs, so derive it from the string instead of trusting the caller's
    context.

    Note this fixes the *value*; displaying it still obeys mp.mp.dps, so use
    `mp.nstr(v, n)` or raise mp.mp.dps to see the digits.
    """
    cleaned = text.strip().replace(" ", "")
    with mp.workdps(_needed_dps(cleaned)):
        return mp.mpf(cleaned)


def _needed_dps(*texts: str) -> int:
    """Decimal precision required to hold these GP-printed numbers exactly."""
    return max(mp.mp.dps, max(sum(c.isdigit() for c in t) for t in texts) + 10)


def _parse_gp_complex(real_text: str, imag_text: str) -> mp.mpc:
    """Build the mpc at full precision, components AND container.

    Parsing the two halves correctly is not enough: `mp.mpc(re, im)` rounds to
    the *global* context, so under the default mp.mp.dps of 15 it truncates
    two perfectly good high-precision mpf values back to doubles.
    """
    real = real_text.strip().replace(" ", "")
    imag = imag_text.strip().replace(" ", "")
    with mp.workdps(_needed_dps(real, imag)):
        return mp.mpc(mp.mpf(real), mp.mpf(imag))


def _reject_unsupported_base(base: GPValue) -> None:
    """Fail fast on bases the Kneser construction here cannot handle.

    Without this, `sexpinit` for a base in (0, 1) never converges: the call sits
    there until `init_timeout` expires, which at the 3600 s default is a silent
    one-hour wait ending in a timeout that names no cause. The README documents
    these bases as unsupported; refusing them up front costs nothing and turns
    an hour into a sentence.

    Only *real* bases are screened. A complex base cannot be classified this
    cheaply, and complex bases in this range are not a documented failure.
    """
    if isinstance(base, str):
        if base.strip() == "e":
            return
        try:
            value = mp.mpf(base.strip())
        except (ValueError, TypeError):
            return          # a GP expression such as "1+I" or "exp(1)"
    elif isinstance(base, (int, float, mp.mpf)):
        value = mp.mpf(base)
    else:
        try:
            z = mp.mpc(base)
        except (ValueError, TypeError):
            return
        if mp.im(z) != 0:
            return
        value = mp.re(z)

    if value <= 0:
        raise ValueError(
            f"base {base!r}: only bases > 0 are supported "
            f"(the construction is undefined for b <= 0).")
    if value < 1:
        raise ValueError(
            f"base {base!r}: bases 0 < b < 1 are not supported -- the "
            f"iteration does not converge and sexpinit would hang until the "
            f"init timeout. Supported: 1 < b < e^(1/e) (regular iteration at "
            f"the attracting fixed point) and b > e^(1/e) (Kneser).")
    if value == 1:
        raise ValueError("base 1 is degenerate (1^x = 1); no tetration exists.")


# Fuehrende/abschliessende Leerzeichen, ein Vorzeichen und die Form ".5"
# kommen in echten Aufrufen vor; ohne sie faellt die Basis stumm in den
# Dezimal-Zweig und reginit lehnt sie erst in GP ab.
_DECIMAL = re.compile(r"^\+?(?:\d+(?:\.\d*)?|\.\d+)$")
# e^(1/e) = 1.44466786100976613365... — enough digits to classify safely
_ETA = Fraction("1.4446678610097661336")


def _exact_if_sub_eta(expr: str) -> str:
    """Send real bases below eta to GP as exact rationals, not decimals.

    The fork's regular-iteration branch (exp-071b) needs the base to roughly
    twice the target precision, because sigma^-1 amplifies every residual by
    lam^-n. A decimal literal is parsed at whatever \\p happens to be and can
    never carry that, so `sexpinit(1.2)` now refuses outright rather than
    returning 15 digits of a 60-digit request. '1.2' and '6/5' are the same
    number, and the rational is exact at any precision.

    Deliberately scoped to the sub-eta interval: every gate base (e, 2, 10,
    1+I, 2+I, 0.8+0.4*I) leaves this function byte-identical, so the Kneser
    path and its measured digit counts cannot shift.

    Applied to STRING bases only. A Python float cannot carry a sub-eta base
    accurately enough -- 1.2 is really 1.1999999999999999555910790149937 --
    and converting its rendered decimal would hand GP the exact rational of
    the wrong number, which reginit would then accept without complaint.
    Floats are left alone so that reginit refuses them audibly.
    """
    stripped = expr.strip()
    if not _DECIMAL.match(stripped):
        return expr
    value = Fraction(stripped)
    if not 1 < value <= _ETA:
        return expr
    return f"{value.numerator}/{value.denominator}"


def _to_gp_number(value: GPValue, digits: int) -> str:
    if isinstance(value, str):
        return value
    # Same trap on the way in: mp.mpc() rounds to the global context, so a
    # high-precision argument passed by a caller who never raised mp.mp.dps
    # would be truncated before it ever reaches GP.
    with mp.workdps(max(mp.mp.dps, digits + 10)):
        z = mp.mpc(value)
        if abs(mp.im(z)) <= mp.mpf(f"1e-{max(20, digits // 2)}"):
            return mp.nstr(mp.re(z), n=digits, min_fixed=0, max_fixed=0)
        real = mp.nstr(mp.re(z), n=digits, min_fixed=0, max_fixed=0)
        imag = mp.nstr(mp.im(z), n=digits, min_fixed=0, max_fixed=0)
        return f"(({real})+({imag})*I)"


@dataclass(frozen=True)
class FatouGPSession:
    gp: "FatouGP"
    base: GPValue

    def eval_batch(self, expressions: Sequence[str]) -> list[mp.mpc]:
        return self.gp.eval_batch(self.base, list(expressions))

    def sexp_batch(self, xs: Sequence[GPValue]) -> list[mp.mpc]:
        return self.gp.sexp_batch(self.base, list(xs))

    def slog_batch(self, ys: Sequence[GPValue]) -> list[mp.mpc]:
        return self.gp.slog_batch(self.base, list(ys))

    def roundtrip_residuals(self, ys: Sequence[GPValue]) -> list[mp.mpc]:
        return self.gp.roundtrip_residuals(self.base, list(ys))

    def sexp(self, x: GPValue) -> mp.mpc:
        return self.gp.sexp(self.base, x)

    def slog(self, y: GPValue) -> mp.mpc:
        return self.gp.slog(self.base, y)


@dataclass
class FatouGP:
    """Persistent PARI/GP tetration engine.

    fatou_gp: engine file; the shortcuts "fork" (optimized fork, ~30x
        faster at high precision — see research/METHODS.md) and
        "original" (unmodified fatou.gp) resolve to the vendored files.
    looplim/nlim: convergence target / iteration cap. Left at None they
        default to looplim=0 (let the engine converge to working
        precision) and nlim=max(30, dps), which is fully converged at
        any precision on BOTH engines — naive calls reproduce the
        published reference values. Expect at least 0.95*dps true
        digits on the fork (measured 302/398/509 at dps 300/400/520); the
        older ~(dps - 24) rule under-promises below dps 500 and
        over-promises above it — see research/METHODS.md, calibration
        section. Only set these explicitly for experiments that
        deliberately under-iterate.
    """
    gp_exe: Path | str | None = None
    fatou_gp: Path | str | None = None
    dps: int = 80
    nlim: int | None = None
    nskip: int = 4
    looplim: int | None = None
    quietmode: int = 1
    persistent: bool = True
    init_timeout: float = 3600.0
    eval_timeout: float = 600.0
    state_cache: bool = True
    n_workers: int = 1
    _workers: dict = field(default_factory=dict, init=False, repr=False, compare=False)
    _pools: dict = field(default_factory=dict, init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        self.gp_exe = find_default_gp_exe() if self.gp_exe is None else Path(self.gp_exe)
        if self.fatou_gp in ("fork", "original"):
            name = "fatou_fork.gp" if self.fatou_gp == "fork" else "fatou.gp"
            self.fatou_gp = Path(__file__).resolve().parent / "vendor" / name
        self.fatou_gp = find_default_fatou_gp() if self.fatou_gp is None else Path(self.fatou_gp)
        if self.looplim is None:
            self.looplim = 0
        if self.nlim is None:
            self.nlim = max(30, self.dps)
        if not self.gp_exe.exists():
            raise FileNotFoundError(f"PARI/GP executable not found: {self.gp_exe}")
        if not self.fatou_gp.exists():
            raise FileNotFoundError(f"fatou.gp not found: {self.fatou_gp}")

    def session(self, base: GPValue) -> FatouGPSession:
        return FatouGPSession(self, base)

    def _base_expr(self, base: GPValue) -> str:
        _reject_unsupported_base(base)
        if isinstance(base, str):
            if base == "e":
                return "exp(1)"
            return _exact_if_sub_eta(base)
        # Bewusst NICHT _exact_if_sub_eta: _to_gp_number() rendert einen
        # float zuerst zu einer Dezimalzeichenkette, und der exakte Bruch
        # DAVON waere der exakte Bruch der falschen Zahl -- Python's 1.2 ist
        # 1.1999999999999999555910790149937... Ein float traegt die Basis
        # prinzipiell nicht genau genug; reginit lehnt sie dann laut ab,
        # statt still die Tetration einer anderen Basis zu liefern.
        return _to_gp_number(base, max(self.dps - 8, 30))

    def _spawn_worker(self, base: GPValue) -> FatouGPWorker:
        if not self.state_cache:
            return self._spawn_full(base)
        key = state_cache.cache_key(
            self._base_expr(base), self.dps, self.nlim, self.nskip,
            self.looplim, self.fatou_gp, self.gp_exe)
        bin_path = state_cache.cache_dir() / f"{key}.gpbin"
        sidecar = state_cache.load_sidecar(bin_path)
        if sidecar is not None:
            worker = self._try_spawn_cached(base, bin_path, sidecar)
            if worker is not None:
                return worker
            state_cache.drop_cache(bin_path)
        worker = self._spawn_full(base)
        try:
            names = state_cache.discover_state_var_names(worker)
            state_cache.dump_state(worker, names, bin_path)
        except (RuntimeError, OSError):
            state_cache.drop_cache(bin_path)  # cache is best-effort
        return worker

    def _spawn_full(self, base: GPValue) -> FatouGPWorker:
        return FatouGPWorker(
            self.gp_exe,
            self._init_lines(base),
            init_timeout=self.init_timeout,
            eval_timeout=self.eval_timeout,
        )

    def _try_spawn_cached(self, base: GPValue, bin_path, sidecar) -> FatouGPWorker | None:
        init = self._init_lines(base)
        # everything except the final sexpinit line, then restore cached state
        lines = init[:-1] + state_cache.restore_lines(sidecar["names"], bin_path)
        try:
            worker = FatouGPWorker(
                self.gp_exe, lines,
                init_timeout=self.init_timeout, eval_timeout=self.eval_timeout)
            anchor = worker.eval([state_cache.ANCHOR_EXPR])[0]
        except (RuntimeError, OSError):
            return None
        expected = mp.mpc(mp.mpf(sidecar["anchor_real"]), mp.mpf(sidecar["anchor_imag"]))
        tolerance = mp.mpf(10) ** (-(min(self.dps, 35)))
        if abs(anchor - expected) > tolerance * max(1, abs(expected)):
            worker.close()
            return None
        return worker

    def _worker_for(self, base: GPValue) -> FatouGPWorker:
        key = self._base_expr(base)
        worker = self._workers.get(key)
        if worker is None or not worker.alive:
            worker = self._spawn_worker(base)
            self._workers[key] = worker
        return worker

    def _pool_for(self, base: GPValue) -> FatouGPPool:
        key = self._base_expr(base)
        pool = self._pools.get(key)
        if pool is None:
            pool = FatouGPPool(lambda: self._spawn_worker(base), self.n_workers)
            self._pools[key] = pool
        return pool

    def close(self) -> None:
        for worker in self._workers.values():
            worker.close()
        self._workers.clear()
        for pool in self._pools.values():
            pool.close()
        self._pools.clear()

    def __enter__(self) -> "FatouGP":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def _init_lines(self, base: GPValue) -> list[str]:
        fatou_path = Path(self.fatou_gp).as_posix()
        return [
            # converged high-dps series overflow the default PARI stack
            "default(parisizemax, 2147483648);",
            f"default(realprecision, {self.dps});",
            f'read("{fatou_path}");',
            f"quietmode={self.quietmode};",
            f"sexpinit({self._base_expr(base)},{self.nlim},{self.nskip},{self.looplim});",
        ]

    def _run_initialized(
        self,
        base: GPValue,
        expressions: list[str],
    ) -> list[mp.mpc]:
        # No `digits` parameter: the callers below have already rendered their
        # arguments with _to_gp_number(x, digits), and the results come back at
        # the engine's working precision. A `digits` argument here was accepted
        # and silently ignored -- results were identical whatever you passed.
        if self.persistent:
            if self.n_workers > 1 and len(expressions) >= 2 * self.n_workers:
                return self._pool_for(base).eval(expressions)
            # an init-phase death (e.g. init_timeout) must NOT be retried
            # blindly — that pays the full init cost twice; only eval-phase
            # deaths get one respawn attempt
            worker = self._worker_for(base)
            try:
                return worker.eval(expressions)
            except WorkerDied:
                self._workers.pop(self._base_expr(base), None)
                return self._worker_for(base).eval(expressions)
        # legacy one-shot path
        lines = self._init_lines(base) + ['print("__BEGIN_RESULTS__")']
        for idx, expr in enumerate(expressions):
            lines.append(f"vv = ({expr});")
            lines.append(f'print("__RES__{idx}")')
            lines.append("print(real(vv));")
            lines.append("print(imag(vv));")
        lines.extend(['print("__END_RESULTS__")', "quit"])
        proc = subprocess.run(
            [str(self.gp_exe), "-q"],
            input="\n".join(lines) + "\n",
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            raise RuntimeError(
                f"PARI/GP failed.\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
            )
        output = proc.stdout.splitlines()
        try:
            start = output.index("__BEGIN_RESULTS__") + 1
            end = output.index("__END_RESULTS__")
        except ValueError as exc:
            raise RuntimeError(f"Could not parse PARI/GP output.\n{proc.stdout}") from exc

        parsed: list[mp.mpc] = []
        cursor = start
        while cursor < end:
            marker = output[cursor].strip()
            if not marker.startswith("__RES__"):
                cursor += 1
                continue
            if cursor + 2 >= end:
                break
            real_text = output[cursor + 1]
            imag_text = output[cursor + 2]
            parsed.append(_parse_gp_complex(real_text, imag_text))
            cursor += 3
        if len(parsed) != len(expressions):
            raise RuntimeError(
                f"Expected {len(expressions)} results, got {len(parsed)}.\nParsed output:\n{proc.stdout}"
            )
        return parsed

    def eval_batch(
        self,
        base: GPValue,
        expressions: list[str],
    ) -> list[mp.mpc]:
        return self._run_initialized(base, expressions)

    def sexp_batch(self, base: GPValue, xs: list[GPValue]) -> list[mp.mpc]:
        digits = max(self.dps - 8, 30)
        expressions = [f"sexp({_to_gp_number(x, digits)})" for x in xs]
        return self._run_initialized(base, expressions)

    def slog_batch(self, base: GPValue, ys: list[GPValue]) -> list[mp.mpc]:
        digits = max(self.dps - 8, 30)
        expressions = [f"slog({_to_gp_number(y, digits)})" for y in ys]
        return self._run_initialized(base, expressions)

    def roundtrip_residuals(self, base: GPValue, ys: list[GPValue]) -> list[mp.mpc]:
        digits = max(self.dps - 8, 30)
        expressions = [
            f"sexp(slog({_to_gp_number(y, digits)})) - ({_to_gp_number(y, digits)})"
            for y in ys
        ]
        return self._run_initialized(base, expressions)

    def sexp(self, base: GPValue, x: GPValue) -> mp.mpc:
        return self.sexp_batch(base, [x])[0]

    def slog(self, base: GPValue, y: GPValue) -> mp.mpc:
        return self.slog_batch(base, [y])[0]
