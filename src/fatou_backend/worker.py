"""Long-lived gp.exe process: init once (sexpinit), then stream evaluations."""
from __future__ import annotations

import queue
import subprocess
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
    return mp.mpf(text.strip().replace(" ", ""))


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
        for stream in (self._proc.stdin, self._proc.stdout):
            try:
                if stream is not None:
                    stream.close()
            except Exception:
                pass

    def __enter__(self) -> "FatouGPWorker":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def __del__(self) -> None:
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
            if "***" in line:
                errors.append(line)

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
            parsed.append(mp.mpc(_parse_gp_scalar(real_text), _parse_gp_scalar(imag_text)))
            cursor += 3
        if len(parsed) != len(expressions):
            self.close()
            raise RuntimeError(
                f"Expected {len(expressions)} results, got {len(parsed)}.\n"
                + "\n".join(output))
        return parsed
