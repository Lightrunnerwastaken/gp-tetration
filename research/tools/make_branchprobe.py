"""Build an instrumented copy of fatou_fork.gp that counts sfunc branch usage.

Measures, per sexpinit run:
  brdirect / brup / brdn  -- how many sfunc samples take the direct / upper
                             theta / lower theta path
  brhorner                -- how many full ct-Horner evaluations (icabel) run
  brshift                 -- how many period shifts the while-loops perform
  brimin                  -- min |imag(y2)| seen (branch-decision margin)
  brtmin                  -- min relative period-tie margin
                             (|y1-y2+-P| - |y1-y2|) / |P|

The margins answer the safety question for F3(i): can the three discrete
decisions that y2 feeds be taken at reduced precision?

No fork mutation -- writes research/tools/fatou_branchprobe.gp.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = Path(__file__).resolve().parent / "fatou_branchprobe.gp"

GLOBALS = """
/* ---- branch probe counters (research/tools/make_branchprobe.py) ---- */
brdirect=0; brup=0; brdn=0; brhorner=0; brshift=0;
brimin=-1; brtmin=-1; brtheta=0;
brreset() = { brdirect=0; brup=0; brdn=0; brhorner=0; brshift=0;
              brimin=-1; brtmin=-1; brtheta=0; }
brreport() = {
  print("BRPROBE direct=", brdirect, " up=", brup, " dn=", brdn,
        " horner=", brhorner, " shift=", brshift);
  print("BRPROBE imin=", brimin, " tmin=", brtmin, " thetacalls=", brtheta);
}
brnote(im, tie) = {
  if ((brimin < 0) || (im < brimin), brimin = im);
  if ((brtmin < 0) || (tie < brtmin), brtmin = tie);
}
"""

# 1) counter in icabel: every call is one O(N) polynomial evaluation
PAT_ICABEL = "icabel(zz) = {\n  local(h, A);"
REP_ICABEL = "icabel(zz) = {\n  local(h, A);\n  brhorner++;"

# 2) direct branch
PAT_DIRECT = "  if ((abs(z-circc)<ircircr)||(thetamode==0),\n    y1 = icabel(z) + n;"
REP_DIRECT = "  if ((abs(z-circc)<ircircr)||(thetamode==0),\n    brdirect++;\n    y1 = icabel(z) + n;"

# 3) upper theta branch: count, and record both decision margins
PAT_UP = """    if (imag(y2)>0,
      if (swon && swidx>0 && swisfv[swidx]==1,"""
REP_UP = """    brtheta++;
    if (imag(y2)>0,
      brup++;
      if (swon && swidx>0 && swisfv[swidx]==1,"""

PAT_UP_W = """      while (abs(y1-y2)>abs(y1-y2+Period), y1=y1+Period);
      while (abs(y1-y2)>abs(y1-y2-Period), y1=y1-Period);
      y1=y1-polcoeff(tht,0);
      z = y1;"""
REP_UP_W = """      while (abs(y1-y2)>abs(y1-y2+Period), y1=y1+Period; brshift++);
      while (abs(y1-y2)>abs(y1-y2-Period), y1=y1-Period; brshift++);
      brnote(abs(imag(y2)),
             min(abs(y1-y2+Period)-abs(y1-y2),
                 abs(y1-y2-Period)-abs(y1-y2))/abs(Period));
      y1=y1-polcoeff(tht,0);
      z = y1;"""

PAT_DN = """    ,
      if (swon && swidx>0 && swisfv[swidx]==2,"""
REP_DN = """    ,
      brdn++;
      if (swon && swidx>0 && swisfv[swidx]==2,"""

PAT_DN_W = """      while (abs(y1-y2)>abs(y1-y2+Period2), y1=y1+Period2);
      while (abs(y1-y2)>abs(y1-y2-Period2), y1=y1-Period2);
      y1=y1-polcoeff(tht2,0);
      z = y1;"""
REP_DN_W = """      while (abs(y1-y2)>abs(y1-y2+Period2), y1=y1+Period2; brshift++);
      while (abs(y1-y2)>abs(y1-y2-Period2), y1=y1-Period2; brshift++);
      brnote(abs(imag(y2)),
             min(abs(y1-y2+Period2)-abs(y1-y2),
                 abs(y1-y2-Period2)-abs(y1-y2))/abs(Period2));
      y1=y1-polcoeff(tht2,0);
      z = y1;"""

PATCHES = [
    (PAT_ICABEL, REP_ICABEL),
    (PAT_DIRECT, REP_DIRECT),
    (PAT_UP, REP_UP),
    (PAT_UP_W, REP_UP_W),
    (PAT_DN, REP_DN),
    (PAT_DN_W, REP_DN_W),
]


def main() -> int:
    text = SRC.read_text(encoding="utf-8", errors="strict")
    for i, (pat, rep) in enumerate(PATCHES):
        n = text.count(pat)
        if n != 1:
            print(f"PATCH {i} matched {n} times -- abort", file=sys.stderr)
            return 1
        text = text.replace(pat, rep, 1)
    text = GLOBALS + text
    DST.write_bytes(text.encode("utf-8"))
    print(f"wrote {DST} ({len(text)} chars)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
