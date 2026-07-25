"""exp-041b: cache the theta branch's y2 per sample for a whole grid stretch.

exp-041a (reduced precision for the same call) gave only 1.05x at dps 200.
This variant removes the call entirely on warm passes, which is the UPPER
BOUND of the whole "y2 is decisions-only" front:

  y2 changes by ~10^-re between passes; the decisions it feeds have measured
  margins 0.795 (branch) and 0.9998 (period tie), and the period loops shift
  0 times. A y2 from the first pass of the stretch therefore yields bitwise
  identical decisions for every later pass of that stretch.

Cache lives in the existing sw* family (allocated and invalidated with
swkey = [samples, w, r]), so a grid change re-derives it.

Writes research/tools/_exp041b.gp.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = Path(__file__).resolve().parent / "_exp041b.gp"

PATCHES: list[tuple[str, str]] = [
    # 1) globals
    ("swon=0; swidx=0; swkey=0; swz=0; swn=0; swvalid=0; swb=0;",
     "swon=0; swidx=0; swkey=0; swz=0; swn=0; swvalid=0; swb=0;\n"
     "/* exp-041b: per-stretch cache of the theta branch's decision-only y2 */\n"
     "swy2=0; swy2v=0;"),
    # 2) allocation alongside the other per-grid caches
    ("""      swisf = vector(samples);
      swisfv = vector(samples);""",
     """      swisf = vector(samples);
      swisfv = vector(samples);
      swy2 = vector(samples);
      swy2v = vector(samples);"""),
    # 3) the call site
    ("    y2 = icabel(zc); /* exp-015: lazy — only the theta path needs it */",
     """    /* exp-015 lazy; exp-041b: y2 only feeds discrete decisions whose
       margins are O(1), so one evaluation per grid stretch suffices. */
    if (swon && swidx>0 && swy2v[swidx],
      y2 = swy2[swidx];
    ,
      y2 = icabel(zc);
      if (swon && swidx>0, swy2[swidx]=y2; swy2v[swidx]=1);
    );"""),
]


def main() -> int:
    text = SRC.read_text(encoding="utf-8")
    for i, (pat, rep) in enumerate(PATCHES):
        if text.count(pat) != 1:
            print(f"patch {i} matched {text.count(pat)} times -- abort", file=sys.stderr)
            return 1
        text = text.replace(pat, rep, 1)
    DST.write_bytes(text.encode("utf-8"))
    print(f"wrote {DST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
