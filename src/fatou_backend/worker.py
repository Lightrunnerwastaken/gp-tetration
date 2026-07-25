"""Long-lived gp.exe process: init once (sexpinit), then stream evaluations."""
from __future__ import annotations

import queue
import subprocess
import sys
import threading
import time
from pathlib import Path

import mpmath as mp


class WorkerDied(RuntimeError):
    """The gp process is gone (crash, EOF, or timeout kill)."""


class _ReaderThread(threading.Thread):
    def __init__(self, stream, out_queue: "queue.Queue[str | None]") -> None:
        super().__init__(daemon=True)
        self._stream = stream
        self._queue = out_queue

    def run(self) -> None:
        for line in iter(self._stream.readline, ""):
            self._queue.put(line.rstrip("\r\n"))
        self._queue.put(None)  # EOF sentinel


def _parse_gp_scalar(text: str) -> mp.mpf:
    """Parse at a precision derived from the string, not from mp.mp.dps.

    See the twin of this function in gp_backend.py: mpmath's global precision
    defaults to 15 dps, so parsing a high-precision GP result under a caller
    who never raised it silently returns a double.
    """
    cleaned = text.strip().replace(" ", "")
    need = max(mp.mp.dps, sum(c.isdigit() for c in cleaned) + 10)
    with mp.workdps(need):
        return mp.mpf(cleaned)


def _parse_gp_complex(real_text: str, imag_text: str) -> mp.mpc:
    """Build the mpc at full precision -- components AND container.

    `mp.mpc(re, im)` rounds to the global context, so parsing the halves
    correctly is not enough; under the default mp.mp.dps of 15 the container
    truncates them straight back to doubles.
    """
    real = real_text.strip().replace(" ", "")
    imag = imag_text.strip().replace(" ", "")
    need = max(mp.mp.dps,
               max(sum(c.isdigit() for c in real),
                   sum(c.isdigit() for c in imag)) + 10)
    with mp.workdps(need):
        return mp.mpc(mp.mpf(real), mp.mpf(imag))


class FatouGPWorker:
    def __init__(self, gp_exe: Path | str, init_lines: list[str],
                 init_timeout: float = 3600.0, eval_timeout: float = 600.0) -> None:
        self.debug_init_lines = list(init_lines)
        self.eval_timeout = eval_timeout
        self._tag = 0
        self._alive = False
        self._proc = subprocess.Popen(
            [str(gp_exe), "-q"],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1,
        )
        self._queue: "queue.Queue[str | None]" = queue.Queue()
        self._reader = _ReaderThread(self._proc.stdout, self._queue)
        self._reader.start()
        self._alive = True
        self._send(init_lines + ['print("__READY__")'])
        self._read_until("__READY__", init_timeout)

    # -- lifecycle ---------------------------------------------------------

    @property
    def alive(self) -> bool:
        return self._alive and self._proc.poll() is None

    def close(self) -> None:
        if getattr(self, "_proc", None) is None:
            return
        self._alive = False
        try:
            if self._proc.poll() is None:
                self._proc.kill()
                self._proc.wait(timeout=10)
        except Exception:
            pass
        try:
            if self._proc.stdin is not None:
                self._proc.stdin.close()
        except Exception:
            pass
        # the reader thread owns the stdout buffer lock while blocked in
        # readline; only close the stream once the thread has seen EOF
        reader = getattr(self, "_reader", None)
        if reader is not None and reader.is_alive() and not sys.is_finalizing():
            reader.join(timeout=5)
        if reader is None or not reader.is_alive():
            try:
                if self._proc.stdout is not None:
                    self._proc.stdout.close()
            except Exception:
                pass

    def __enter__(self) -> "FatouGPWorker":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def __del__(self) -> None:
        if sys.is_finalizing():
            return  # interpreter teardown: daemon threads are frozen, OS reclaims
        self.close()

    # -- protocol ----------------------------------------------------------

    def _send(self, lines: list[str]) -> None:
        if not self.alive:
            raise WorkerDied("gp process is not running")
        try:
            self._proc.stdin.write("\n".join(lines) + "\n")
            self._proc.stdin.flush()
        except OSError as exc:
            self.close()
            raise WorkerDied(f"failed to write to gp: {exc}") from exc

    def _read_until(self, sentinel: str, timeout: float) -> list[str]:
        deadline = time.perf_counter() + timeout
        collected: list[str] = []
        errors: list[str] = []
        while True:
            remaining = deadline - time.perf_counter()
            if remaining <= 0:
                self.close()
                raise WorkerDied(
                    f"timeout after {timeout}s waiting for {sentinel!r}. "
                    f"Last output:\n" + "\n".join(collected[-20:]))
            try:
                line = self._queue.get(timeout=min(remaining, 1.0))
            except queue.Empty:
                continue
            if line is None:
                self.close()
                raise WorkerDied(
                    "gp process ended unexpectedly. Last output:\n"
                    + "\n".join(collected[-20:]))
            if line.strip() == sentinel:
                if errors:
                    self.close()
                    raise RuntimeError("PARI/GP error:\n" + "\n".join(errors))
                return collected
            collected.append(line)
            if "***" in line and "Warning" not in line:
                errors.append(line)  # PARI warnings (e.g. parisizemax) are not errors

    def eval_raw(self, lines: list[str], sentinel: str, timeout: float) -> list[str]:
        self._send(lines + [f'print("{sentinel}")'])
        return self._read_until(sentinel, timeout)

    def eval(self, expressions: list[str]) -> list[mp.mpc]:
        self._tag += 1
        sentinel = f"__END__{self._tag}"
        lines: list[str] = []
        for idx, expr in enumerate(expressions):
            lines.append(f"vv = ({expr});")
            lines.append(f'print("__RES__{idx}")')
            lines.append("print(real(vv));")
            lines.append("print(imag(vv));")
        output = self.eval_raw(lines, sentinel, self.eval_timeout)

        parsed: list[mp.mpc] = []
        cursor = 0
        while cursor < len(output):
            marker = output[cursor].strip()
            if not marker.startswith("__RES__"):
                cursor += 1
                continue
            if cursor + 2 > len(output) - 1:
                break
            real_text = output[cursor + 1]
            imag_text = output[cursor + 2]
            parsed.append(_parse_gp_complex(real_text, imag_text))
            cursor += 3
        if len(parsed) != len(expressions):
            self.close()
            raise RuntimeError(
                f"Expected {len(expressions)} results, got {len(parsed)}.\n"
                + "\n".join(output))
        return parsed
