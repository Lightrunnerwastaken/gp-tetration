"""Fixed benchmark/gate workload. Frozen once reference values exist."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Case:
    kind: str   # "sexp" | "slog"
    base: str   # GP base expression, e.g. "e", "2", "1+I"
    arg: str    # GP-parsable argument literal
    dps: int
    tier: str   # "throughput" | "precision" | "anchor"


THROUGHPUT_BASES = ("e", "2", "10")
THROUGHPUT_DPS = 80
# 16 sexp arguments in [-0.4375, 1.4375]
SEXP_ARGS = tuple(str(-7 / 16 + j / 8) for j in range(16))
# 16 slog arguments in [1.25, 5.0]
SLOG_ARGS = tuple(str(1.25 + j / 4) for j in range(16))
PRECISION_CASES = (("e", "0.5"), ("2", "0.5"))
PRECISION_DPS = (200,)
PRECISION_DPS_DEEP = (200, 500, 1000)
ANCHOR_BASES = ("1+I", "0.8+0.4*I", "2+I")
ANCHOR_DPS = 80


def all_cases(deep: bool = False) -> list[Case]:
    cases: list[Case] = []
    for base in THROUGHPUT_BASES:
        for arg in SEXP_ARGS:
            cases.append(Case("sexp", base, arg, THROUGHPUT_DPS, "throughput"))
        for arg in SLOG_ARGS:
            cases.append(Case("slog", base, arg, THROUGHPUT_DPS, "throughput"))
    for base, arg in PRECISION_CASES:
        for dps in (PRECISION_DPS_DEEP if deep else PRECISION_DPS):
            cases.append(Case("sexp", base, arg, dps, "precision"))
    for base in ANCHOR_BASES:
        cases.append(Case("sexp", base, "0.5", ANCHOR_DPS, "anchor"))
    return cases


def case_id(c: Case) -> str:
    return f"{c.kind}|{c.base}|{c.arg}|{c.dps}"
