"""Build a section-profiled copy of the CURRENT fatou_fork.gp.

An older profiled copy exists (research/tools/fatou_profile6.gp) but it was
patched from a pre-exp-032/037 fork, so its numbers no longer describe the
engine that runs today (exp-040 lesson: re-measure before every mutation).

Per loop() iteration it prints one PROF line with:

    tth   ms in thtaylor (theta rebuild, both maps)
    tsmp  ms in staylor's sampling loop  (the N x Horner block)
    text  ms in staylor's extraction     (FFT / Bluestein)
    tsta  ms in staylor overall
    icfull / exinc / thfull   cold-vs-warm flags for ct / extraction / theta

Cold passes (icfull=1) are grid-stretch starts: full Horner at full working
precision. Warm passes evaluate only the diff at icdig digits. The split
between those two is the number that decides where an optimization must hit.

Writes research/tools/_prof.gp. No fork mutation.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = Path(__file__).resolve().parent / "_prof.gp"

PATCHES: list[tuple[str, str]] = [
    # globals
    ("quietmode=0;\n/* I added || (real(Period)>47)",
     "quietmode=0;\npfsmp=0; pfext=0; pfa=0; pfb=0; pftth=0; pfsta=0; pft0=0; pfs0=0;\n"
     "/* I added || (real(Period)>47)"),
    # staylor: sampling loop
    ("""  if (complextaylor,
    for(s=1, samples,
      x1=-1+-1/(samples)+(2*s/samples);""",
     """  pfa=getabstime();
  if (complextaylor,
    for(s=1, samples,
      x1=-1+-1/(samples)+(2*s/samples);"""),
    ("""  swon = 0;
  swidx = 0;
  wtaylor=0;""",
     """  pfsmp=getabstime()-pfa;
  swon = 0;
  swidx = 0;
  wtaylor=0;"""),
    # staylor: extraction block
    ("""  if (subeta==0,
    /* exp-027: on an unchanged grid, extract only the diff of the samples""",
     """  pfb=getabstime();
  if (subeta==0,
    /* exp-027: on an unchanged grid, extract only the diff of the samples"""),
    ("""  if (st==0, st=terms);
  wtaylor=precision(wtaylor,precis);""",
     """  pfext=getabstime()-pfb;
  if (st==0, st=terms);
  wtaylor=precision(wtaylor,precis);"""),
    # loop(): theta rebuild timer
    ("""    if (thetamode,
      thetamode++;  /* used by renormr */
      thsamples=floor(re*tht_re_mult+6);""",
     """    pft0=getabstime();
    if (thetamode,
      thetamode++;  /* used by renormr */
      thsamples=floor(re*tht_re_mult+6);"""),
    ("""    ctsamples = floor((stopterms+20)*1.03);""",
     """    pftth=getabstime()-pft0;
    ctsamples = floor((stopterms+20)*1.03);"""),
    # loop(): staylor timer + PROF line
    ("""    z=staylor(circc,r,ctsamples); ct=z[1]; stopterms=z[2];
    rr=renormr(ct);
    relast=re;""",
     """    pfs0=getabstime();
    z=staylor(circc,r,ctsamples); ct=z[1]; stopterms=z[2];
    pfsta=getabstime()-pfs0;
    print("PROF n=", n, " re=", precision(re,6), " cts=", ctsamples,
          " ths=", thsamples, " dps=", default(realprecision),
          " icdig=", icdig, " tth=", pftth, " tsmp=", pfsmp,
          " text=", pfext, " tsta=", pfsta,
          " icfull=", icfull, " exinc=", exinc, " thfull=", thfull);
    rr=renormr(ct);
    relast=re;"""),
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
