"""exp-051 candidate: the theta grid is unquantized below 64, so its caches
never engage in a large part of the run.

exp-026 quantizes the theta grid to 64-steps so it stays put for many
iterations -- that is the precondition for its per-index superf cache and for
the incremental Horner (thdct at reduced precision). But the guard is

    ... && (samples > 64), samples = 64*ceil(samples/64)

so BELOW 64 the grid changes on essentially every iteration, thfull stays 1,
the caches are reallocated and everything is recomputed at full precision.

Measured (base e, dps 300, research/tools/make_profile.py): thfull = 1 on 60 of
232 iterations, and they are exactly the contiguous block n = 37..108 where
thsamples runs 27..62 -- the entire sub-64 range. That block is not cheap: the
ct degree there is already several hundred.

Fix: quantize below 64 as well, with a smaller step so the overshoot stays
small. thsamples grows ~0.5 per iteration, so an 8-step grid is stable for
~16 iterations -- enough for the incremental path to pay off.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
HERE = Path(__file__).resolve().parent

PAT = ("  if ((subeta==0) && (efam || (complextaylor==0)) && (n==1) && (samples>64), "
       "samples = 64*ceil(samples/64));")
REP = """  /* exp-051: exp-026 quantized only above 64, so below it the theta grid
     moved every iteration and thfull stayed 1 -- caches reallocated and the
     incremental path never engaged. Measured at dps 300: that was 60 of 232
     iterations, the contiguous block thsamples 27..62. Quantize there too,
     with a smaller step (%(Q)d) so the overshoot stays small. */
  if ((subeta==0) && (efam || (complextaylor==0)) && (n==1),
    if (samples > 64,
      samples = 64*ceil(samples/64)
    ,
      if (samples > %(Q)d, samples = %(Q)d*ceil(samples/%(Q)d))));"""


def build(q: int, dst: Path) -> Path:
    text = SRC.read_text(encoding="utf-8")
    if text.count(PAT) != 1:
        print(f"pattern matched {text.count(PAT)} times -- abort", file=sys.stderr)
        raise SystemExit(1)
    dst.write_bytes(text.replace(PAT, REP % {"Q": q}, 1).encode("utf-8"))
    return dst


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quantum", type=int, default=8)
    a = ap.parse_args()
    print(f"wrote {build(a.quantum, HERE / f'_exp051_{a.quantum}.gp')} (sub-64 quantum {a.quantum})")
