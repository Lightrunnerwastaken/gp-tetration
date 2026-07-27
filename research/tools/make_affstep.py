"""Build an affstep harness from the CURRENT fork.

research/tools/fatou_e3diag.gp carries the same machinery but is a fork copy
from 2026-07-12 -- nine keeps old, with ctr = 0.72 where the shipped engine has
81/125 = 0.648. Every structural measurement taken through it describes an
engine that no longer exists, which is exactly the trap the 2026-07-27 audit
flagged. This patcher regenerates the harness from whatever fork is shipped
today, so the numbers describe the artifact.

What it adds:
  affstop      stop the Picard loop as soon as re >= affstop (freeze a state)
  affcts/affths  the grid sizes at the moment of the freeze
  affstep()    ONE loop iteration as an operator on the global ct, at those
               FIXED grids -- the branch-locked step G, so that
               G(c + d) - G(c) = A*d exactly (the iteration is affine on a
               fixed grid, verified to 4.5e-81 in E3).

Note affstep() is only bit-exactly deterministic from the SECOND call onward
(measured: 3rd minus 2nd = exactly 0, 2nd minus 1st = 5e-219); the incremental
caches warm on the first. Any probe must discard one call before taking its
reference.

Writes research/tools/_affstep.gp.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = Path(__file__).resolve().parent / "_affstep.gp"

PATCHES: list[tuple[str, str]] = [
    ("quietmode=0;\n/* I added || (real(Period)>47)",
     "quietmode=0;\naffstop=0; affcts=0; affths=0;\n"
     "/* I added || (real(Period)>47)"),
    # freeze point: right after re is updated for this iteration
    ("""    if ((n>3) && ((re-relast) < 0.1), nskip--);
    if (quietmode==0,""",
     """    if ((n>3) && ((re-relast) < 0.1), nskip--);
    affcts = ctsamples; affths = thsamples;
    if ((affstop > 0) && (re >= affstop), break());
    if (quietmode==0,"""),
    # the operator itself, placed just before sexpinit
    ("""sexpinit(b,nlim,nskip,looplim) = {""",
     """/* one loop iteration as an operator on the global ct, at the frozen grids */
affstep() = {
  local(z);
  ct = precision(ct, precis);
  ct = ct + renormr(ct);
  if (thetamode,
    thetamode++;
    tht = thtaylor(1, affths);
    if (complextaylor, tht2 = thtaylor(2, affths));
  );
  z = staylor(circc, if (thetamode, circr*ctr, circr), affcts);
  ct = z[1];
  return(ct);
}

sexpinit(b,nlim,nskip,looplim) = {"""),
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
