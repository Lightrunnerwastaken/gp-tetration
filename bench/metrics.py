"""Digit-accuracy metric and reference loading. Pure functions, no GP calls."""
from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp


def correct_digits(value: mp.mpc, reference: mp.mpc, cap: int) -> float:
    err = abs(mp.mpc(value) - mp.mpc(reference))
    if err == 0:
        return float(cap)
    scale = max(mp.mpf(1), abs(mp.mpc(reference)))
    digits = -mp.log10(err / scale)
    return float(min(mp.mpf(cap), digits))


def load_reference(path: Path) -> dict[str, mp.mpc]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    return {
        key: mp.mpc(mp.mpf(entry["real"]), mp.mpf(entry["imag"]))
        for key, entry in payload["values"].items()
    }
