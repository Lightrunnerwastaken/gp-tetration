"""Count mixdft twiddle-cache hits vs misses. No timing, exact counts.

The exp-061 cache keys on (nn, realprecision) and accepts a slot only when it
was built at a precision >= the current one:

    if ((mx1n == nn) && (mx1p >= default(realprecision)), sl = 1);

But loop() raises the global working precision on EVERY iteration
(`default(realprecision, max(48, min(precis, 2*floor(re)+60)))`, line ~1737),
so while the ladder is still climbing, a slot built at iteration i is stale at
iteration i+1 by construction. Only once the ladder saturates at `precis` does
the precision stop rising -- and inside the extraction block it then *falls*
(exdig = realprecision - re + 40), at which point every lookup should hit.

So the cache is expected to be useless in the first half of a run and free
money in the second. That is a guess about a data structure that costs ~2*nn
full-precision numbers of memory, which at dps 1020 / nn 7168 is not nothing.
Measure it instead: this prints one line per mixdft call.

Writes research/tools/_mxcount.gp. No fork mutation.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = Path(__file__).resolve().parent / "_exp061.gp"
DST = Path(__file__).resolve().parent / "_mxcount.gp"

PATCHES: list[tuple[str, str]] = [
    ("mx1n=0; mx1p=0; mx1t=0; mx1w=0; mx1o=0; mx2n=0; mx2p=0; mx2t=0; mx2w=0; mx2o=0;",
     "mx1n=0; mx1p=0; mx1t=0; mx1w=0; mx1o=0; mx2n=0; mx2p=0; mx2t=0; mx2w=0; mx2o=0;\n"
     "mxhit=0; mxmiss=0;"),
    ("""  if (sl == 0,
    w = geopow(exp(-2*Pi*I/m), m-1);""",
     """  if (sl, mxhit++, mxmiss++);
  print("MX nn=", nn, " prec=", default(realprecision),
        " slot=", sl, " hit=", mxhit, " miss=", mxmiss);
  if (sl == 0,
    w = geopow(exp(-2*Pi*I/m), m-1);"""),
]


def main() -> int:
    if not SRC.exists():
        print(f"{SRC} missing -- run make_exp061.py first", file=sys.stderr)
        return 1
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
