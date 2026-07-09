"""Round-robin pool of persistent workers for large batches."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from typing import Callable

import mpmath as mp

from .worker import FatouGPWorker


class FatouGPPool:
    def __init__(self, spawn: Callable[[], FatouGPWorker], n_workers: int) -> None:
        if n_workers < 1:
            raise ValueError("n_workers must be >= 1")
        self._workers = [spawn() for _ in range(n_workers)]
        self._executor = ThreadPoolExecutor(max_workers=n_workers)

    def eval(self, expressions: list[str]) -> list[mp.mpc]:
        n = len(self._workers)
        chunks = [expressions[i::n] for i in range(n)]
        futures = [
            self._executor.submit(worker.eval, chunk)
            for worker, chunk in zip(self._workers, chunks) if chunk
        ]
        chunk_results = [f.result() for f in futures]
        merged: list[mp.mpc | None] = [None] * len(expressions)
        for i, results in enumerate(chunk_results):
            for j, value in enumerate(results):
                merged[i + j * n] = value
        return merged  # type: ignore[return-value]

    def close(self) -> None:
        for worker in self._workers:
            worker.close()
        self._executor.shutdown(wait=False)
