"""Build a PHASE-SPLIT profiler for thtaylor/thfunc.

thtaylor's sampling loop has no cross-sample dependencies (every cache slot is
indexed by s), so the loop can be split into independent passes without
changing a single arithmetic operation.  Each pass is timed separately, which
gives an exact decomposition of the theta block:

  thp  prep      thdct = ct-icct, vecmax, precision()  (once per call)
  th1  grid      x1 -> tcrc[s]=exp(Pi*I*x1) -> y=log(tcrc)/(2 Pi I)
  th2  superf    superf(zth+y) or the exp-026 cache hit
  th3  argred    precision(p-circc, thdig)
  th4  horner    subst(ct|thdct, x, .)
  th5  branchA   the two logs of abelest's closed form (cold only)
  th6  assemble  (A+h)-y-zth and the cache stores
  th7  extract   fft/mixdft + powers + Polrev + precision

Writes research/tools/_thphase.gp.  No fork mutation.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = Path(__file__).resolve().parent / "_thphase.gp"

OLD_LOOP = """  for(s=1, samples,
    x1 = -1 + -1/(samples) + (2*s/samples); /* -Pi to Pi */
    tcrc[s] = exp(Pi*I*x1);
    thidx = s;
    t_est[s]= thfunc(tcrc[s],n);
  );
  thon = 0; thidx = 0;
"""

NEW_LOOP = """  if (thon,
    thyv = vector(samples); thpv = vector(samples);
    thqv = vector(samples); thhv = vector(samples);
    pfc=getabstime();
    for(s=1, samples,
      x1 = -1 + -1/(samples) + (2*s/samples);
      tcrc[s] = exp(Pi*I*x1);
      thyv[s] = log(tcrc[s])/(2*Pi*I);
    );
    pf1=pf1+(getabstime()-pfc);
    pfc=getabstime();
    for(s=1, samples,
      if (thswv[s],
        thpv[s] = thsw[s];
      ,
        thpv[s] = superf(zth+thyv[s]);
        thsw[s] = thpv[s]; thswv[s] = 1;
      );
    );
    pf2=pf2+(getabstime()-pfc);
    pfc=getabstime();
    if (thfull,
      for(s=1, samples, thqv[s] = thpv[s]-circc);
    ,
      for(s=1, samples, thqv[s] = precision(thpv[s]-circc, thdig));
    );
    pf3=pf3+(getabstime()-pfc);
    pfc=getabstime();
    if (thfull,
      for(s=1, samples, thhv[s] = subst(ct, x, thqv[s]));
    ,
      for(s=1, samples, thhv[s] = thicv[s] + subst(thdct, x, thqv[s]));
    );
    pf4=pf4+(getabstime()-pfc);
    pfc=getabstime();
    if (thfull,
      for(s=1, samples,
        thicA[s] = rlnlm*(log(I*(thpv[s]-L))-Pi*I/2)
                 + rlnlm2*(log(-I*(thpv[s]-L2))+Pi*I/2) + sfunczero);
    );
    pf5=pf5+(getabstime()-pfc);
    pfc=getabstime();
    for(s=1, samples,
      thicv[s] = thhv[s];
      t_est[s] = (thicA[s] + thhv[s]) - thyv[s] - zth;
    );
    pf6=pf6+(getabstime()-pfc);
  ,
    for(s=1, samples,
      x1 = -1 + -1/(samples) + (2*s/samples); /* -Pi to Pi */
      tcrc[s] = exp(Pi*I*x1);
      thidx = s;
      t_est[s]= thfunc(tcrc[s],n);
    );
  );
  thon = 0; thidx = 0;
  pfc=getabstime();
"""

PATCHES: list[tuple[str, str]] = [
    ("quietmode=0;\n/* I added || (real(Period)>47)",
     "quietmode=0;\npf1=0;pf2=0;pf3=0;pf4=0;pf5=0;pf6=0;pf7=0;pfp=0;pfc=0;pfsta=0;pfa=0;\n"
     "thyv=0;thpv=0;thqv=0;thhv=0;\n"
     "/* I added || (real(Period)>47)"),
    # prep block of thtaylor
    ("""  wtaylor=0;
  if ((subeta==0) && (efam || (complextaylor==0)) && (n==1),
    thfull = 1;""",
     """  wtaylor=0;
  pfc=getabstime();
  if ((subeta==0) && (efam || (complextaylor==0)) && (n==1),
    thfull = 1;"""),
    ("""    thskey = samples;
    thon = 1;
  );
""",
     """    thskey = samples;
    thon = 1;
  );
  pfp=pfp+(getabstime()-pfc);
"""),
    (OLD_LOOP, NEW_LOOP),
    ("""  wtaylor=Polrev(cf);
  wtaylor=precision(wtaylor,precis);
  return(wtaylor);
}""",
     """  wtaylor=Polrev(cf);
  wtaylor=precision(wtaylor,precis);
  pf7=pf7+(getabstime()-pfc);
  return(wtaylor);
}"""),
    # staylor total, for the denominator
    ("""staylor( w,r,samples) = {
  local(rinv,""",
     """staylor( w,r,samples) = {
  pfa=getabstime();
  local(rinv,"""),
    ("""  if (st==0, st=terms);
  wtaylor=precision(wtaylor,precis);
  if (complextaylor, wtaylor=subst(wtaylor,x,x*exp(-argc*I)));
  return([wtaylor,st]);""",
     """  if (st==0, st=terms);
  wtaylor=precision(wtaylor,precis);
  if (complextaylor, wtaylor=subst(wtaylor,x,x*exp(-argc*I)));
  pfsta=pfsta+(getabstime()-pfa);
  return([wtaylor,st]);"""),
]


def main() -> None:
    text = SRC.read_text(encoding="utf-8", errors="replace")
    for old, new in PATCHES:
        if old not in text:
            sys.exit(f"PATCH NOT FOUND:\n{old[:250]}")
        text = text.replace(old, new, 1)
    DST.write_text(text, encoding="utf-8")
    print(f"wrote {DST}")


if __name__ == "__main__":
    main()
