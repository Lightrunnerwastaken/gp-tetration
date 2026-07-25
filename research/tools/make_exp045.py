"""exp-045 candidate: re-tune the e-family sampling radius.

exp-023 (2026-07-11) set `if (efam, ctr = ctr*9/10)` after a scan at dps
150-400, and a rescan on 2026-07-12 confirmed it. Since then the cost
structure changed underneath that scan: exp-024 (isuperf cache), exp-026/027
(theta quantization + incremental extraction), exp-030..034 (chirp caches,
unified extraction), exp-037 (512-quantum), exp-041b (theta branch no longer
evaluates ct) and exp-042b (theta extraction via FFT). The dominant N^2
sampling term now carries a larger share, so the product optimum moves to a
SMALLER radius.

Re-scan with the existing ctrmul knob (research/tools/ctr_scan.py, no
mutation), base e:

    dps 200:  ctrmul 1 -> 20.00 s,  0.9 -> 18.99 s,  0.8 -> 18.80 s
    dps 300:  ctrmul 1 -> 72.41 s,  0.9 -> 64.70 s,  0.8 -> 61.34 s
    anchors agree to 203-297 digits, i.e. accuracy is untouched

The gain grows with depth (1.064x at dps 200, 1.181x at dps 300) because the
grid shrinks quadratically while the iteration count only grows linearly.

This mutation folds the scanned factor into the shipped radius:
9/10 -> 9/10 * 4/5 = 18/25. Only the e-family is touched; the fast bases keep
their delicate sample/terms co-evolution, exactly as exp-023 intended.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
HERE = Path(__file__).resolve().parent

PAT = "  if (efam, ctr = ctr*9/10);"
REP = """  /* exp-045: re-scan of the exp-023 radius after ~10 keeps changed the
     cost balance (exp-024/026/027/030-034/037/041b/042b). The N^2 sampling
     term now dominates enough that a smaller radius wins despite the slower
     rate: measured 1.064x at dps 200 and 1.181x at dps 300 with unchanged
     digits, and the gain grows with depth. 9/10 -> %(f)s. */
  if (efam, ctr = ctr*%(f)s);"""


def build(factor: str, dst: Path) -> Path:
    text = SRC.read_text(encoding="utf-8")
    if text.count(PAT) != 1:
        print(f"pattern matched {text.count(PAT)} times -- abort", file=sys.stderr)
        raise SystemExit(1)
    dst.write_bytes(text.replace(PAT, REP % {"f": factor}, 1).encode("utf-8"))
    return dst


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--factor", default="18/25", help="replaces 9/10 for efam")
    a = ap.parse_args()
    tag = a.factor.replace("/", "_")
    print(f"wrote {build(a.factor, HERE / f'_exp045_{tag}.gp')} (efam ctr factor {a.factor})")
