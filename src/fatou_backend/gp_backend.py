from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import mpmath as mp


GPValue = mp.mpf | mp.mpc | float | complex | int | str


def _first_existing(candidates: Iterable[Path]) -> Path | None:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def find_default_gp_exe() -> Path:
    env = os.getenv("FATOU_GP_EXE")
    if env:
        path = Path(env)
        if path.exists():
            return path
    candidates = (
        Path(r"C:\Program Files (x86)\Pari32-2-17-3\gp.exe"),
        Path(r"C:\Program Files\Pari64-2-17-3\gp.exe"),
    )
    path = _first_existing(candidates)
    if path is None:
        raise FileNotFoundError("Could not locate gp.exe. Set FATOU_GP_EXE.")
    return path


def find_default_fatou_gp() -> Path:
    env = os.getenv("FATOU_GP_FILE")
    if env:
        path = Path(env)
        if path.exists():
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
    cleaned = text.strip().replace(" ", "")
    return mp.mpf(cleaned)


def _to_gp_number(value: GPValue, digits: int) -> str:
    if isinstance(value, str):
        return value
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

    def eval_batch(self, expressions: Sequence[str], digits: int | None = None) -> list[mp.mpc]:
        return self.gp.eval_batch(self.base, list(expressions), digits=digits)

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
    gp_exe: Path | str | None = None
    fatou_gp: Path | str | None = None
    dps: int = 80
    nlim: int = 30
    nskip: int = 4
    looplim: int = 35
    quietmode: int = 1

    def __post_init__(self) -> None:
        self.gp_exe = find_default_gp_exe() if self.gp_exe is None else Path(self.gp_exe)
        self.fatou_gp = find_default_fatou_gp() if self.fatou_gp is None else Path(self.fatou_gp)
        if not self.gp_exe.exists():
            raise FileNotFoundError(f"PARI/GP executable not found: {self.gp_exe}")
        if not self.fatou_gp.exists():
            raise FileNotFoundError(f"fatou.gp not found: {self.fatou_gp}")

    def session(self, base: GPValue) -> FatouGPSession:
        return FatouGPSession(self, base)

    def _base_expr(self, base: GPValue) -> str:
        if isinstance(base, str):
            if base == "e":
                return "exp(1)"
            return base
        return _to_gp_number(base, max(self.dps - 8, 30))

    def _run_initialized(
        self,
        base: GPValue,
        expressions: list[str],
        digits: int | None = None,
    ) -> list[mp.mpc]:
        if digits is None:
            digits = max(self.dps - 8, 30)
        fatou_path = self.fatou_gp.as_posix()
        lines = [
            f"\\\\p {self.dps}",
            f'read("{fatou_path}")',
            f"quietmode={self.quietmode};",
            f"sexpinit({self._base_expr(base)},{self.nlim},{self.nskip},{self.looplim});",
            'print("__BEGIN_RESULTS__")',
        ]
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
            parsed.append(mp.mpc(_parse_gp_scalar(real_text), _parse_gp_scalar(imag_text)))
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
        digits: int | None = None,
    ) -> list[mp.mpc]:
        return self._run_initialized(base, expressions, digits=digits)

    def sexp_batch(self, base: GPValue, xs: list[GPValue]) -> list[mp.mpc]:
        digits = max(self.dps - 8, 30)
        expressions = [f"sexp({_to_gp_number(x, digits)})" for x in xs]
        return self._run_initialized(base, expressions, digits=digits)

    def slog_batch(self, base: GPValue, ys: list[GPValue]) -> list[mp.mpc]:
        digits = max(self.dps - 8, 30)
        expressions = [f"slog({_to_gp_number(y, digits)})" for y in ys]
        return self._run_initialized(base, expressions, digits=digits)

    def roundtrip_residuals(self, base: GPValue, ys: list[GPValue]) -> list[mp.mpc]:
        digits = max(self.dps - 8, 30)
        expressions = [
            f"sexp(slog({_to_gp_number(y, digits)})) - ({_to_gp_number(y, digits)})"
            for y in ys
        ]
        return self._run_initialized(base, expressions, digits=digits)

    def sexp(self, base: GPValue, x: GPValue) -> mp.mpc:
        return self.sexp_batch(base, [x])[0]

    def slog(self, base: GPValue, y: GPValue) -> mp.mpc:
        return self.slog_batch(base, [y])[0]
