"""Build a COUNTING probe of the theta path (no timing, exact counts).

Answers: what does the thtaylor/thfunc/superf/isuperf block actually DO per
call -- how many superf calls, how long the fs()-ladder inside superf is, how
many Horner terms at which precision, how many cold vs warm passes.

Writes research/tools/_thcount.gp. No fork mutation.
"""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
DST = Path(__file__).resolve().parent / "_thcount.gp"

PATCHES: list[tuple[str, str]] = [
    ("quietmode=0;\n/* I added || (real(Period)>47)",
     "quietmode=0;\n"
     "cthc=0; ctsupf=0; ctsupn=0; ctsupd=0; ctfull=0; ctwarm=0;\n"
     "cthdeg=0; cthdig=0; ctisf=0; ctisfn=0; ctisfd=0;\n"
     "ctsmp=0; ctsdeg=0; ctsdig=0; ctlvl=0; ctit=0;\n"
     "ctrad=0; ctradn=0; ctradmin=99; ctradmax=0;\n"
     "ctsupmode=0;\n"
     "/* I added || (real(Period)>47)"),

    # ---- superf: count ladder length and series length -------------------
    ("""superf(z) = {
  local(n,i,y);
  if (repelling,
    y=lambda1^(z-superfk);
    n=superfk;
    while (abs(y)>superfr, y=y/(lambda1);n++);
    y = subst(fisl,x,y);
    for (i=1,n,y=fs(y));
    return(y);""",
     """superf(z) = {
  local(n,i,y);
  if (repelling,
    y=lambda1^(z-superfk);
    n=superfk;
    while (abs(y)>superfr, y=y/(lambda1);n++);
    if (ctsupmode, ctsupn=ctsupn+n; ctsupd=ctsupd+length(fisl));
    y = subst(fisl,x,y);
    for (i=1,n,y=fs(y));
    return(y);"""),

    # ---- isuperf: count ladder length ------------------------------------
    ("""isuperf(z) = {
  local(n,y,y1);
  y=z;
  n=0;
  if (repelling,
    while (abs(y-L)>isuperfr,""",
     """isuperf(z) = {
  local(n,y,y1);
  ctisf++;
  y=z;
  n=0;
  if (repelling,
    while (abs(y-L)>isuperfr,"""),
    ("""  y = subst(fsl,x,y-L);
  y = log(y)*rlnlm + n;
  while (abs(imag(y))>abs(imag(y+Period)),y=y+Period);""",
     """  ctisfn=ctisfn+abs(n); ctisfd=ctisfd+length(fsl);
  y = subst(fsl,x,y-L);
  y = log(y)*rlnlm + n;
  while (abs(imag(y))>abs(imag(y+Period)),y=y+Period);"""),

    # ---- thfunc: count calls / superf / Horner degree*precision ----------
    ("""      if (thswv[thidx],
        p = thsw[thidx];
      ,
        p = superf(zth+y);
        thsw[thidx] = p; thswv[thidx] = 1;
      );
      if (thfull,
        A = rlnlm*(log(I*(p-L))-Pi*I/2) + rlnlm2*(log(-I*(p-L2))+Pi*I/2) + sfunczero;
        h = subst(ct, x, (p-circc));
        thicA[thidx] = A; thicv[thidx] = h;
      ,
        h = thicv[thidx] + subst(thdct, x, precision(p-circc, thdig));
        thicv[thidx] = h;
        A = thicA[thidx];
      );""",
     """      cthc++;
      if (thswv[thidx],
        p = thsw[thidx];
      ,
        ctsupf++; ctsupmode=1;
        p = superf(zth+y);
        ctsupmode=0;
        thsw[thidx] = p; thswv[thidx] = 1;
      );
      ctrad=ctrad+abs(p-circc); ctradn++;
      if (abs(p-circc)<ctradmin, ctradmin=abs(p-circc));
      if (abs(p-circc)>ctradmax, ctradmax=abs(p-circc));
      if (thfull,
        ctfull++;
        cthdeg=cthdeg+length(ct); cthdig=cthdig+length(ct)*default(realprecision);
        A = rlnlm*(log(I*(p-L))-Pi*I/2) + rlnlm2*(log(-I*(p-L2))+Pi*I/2) + sfunczero;
        h = subst(ct, x, (p-circc));
        thicA[thidx] = A; thicv[thidx] = h;
      ,
        ctwarm++;
        cthdeg=cthdeg+length(thdct); cthdig=cthdig+length(thdct)*thdig;
        h = thicv[thidx] + subst(thdct, x, precision(p-circc, thdig));
        thicv[thidx] = h;
        A = thicA[thidx];
      );"""),

    # ---- staylor sampling: count Horner work in icabel -------------------
    ("""    if (icfull,
      A = rlnlm*(log(I*(zz-L))-Pi*I/2) + rlnlm2*(log(-I*(zz-L2))+Pi*I/2) + sfunczero;
      h = subst(if (lv, icdl[lv], ct), x, (zz-circc));
      icA[swidx] = A;
      icvals[swidx] = h;
    ,
      h = icvals[swidx] + subst(if (lv, icdl[lv], icdct), x, precision(zz-circc, icdig));
      icvals[swidx] = h;
      A = icA[swidx];
    );""",
     """    if (icfull,
      ctsmp++;
      ctsdeg=ctsdeg+length(if (lv, icdl[lv], ct));
      ctsdig=ctsdig+length(if (lv, icdl[lv], ct))*default(realprecision);
      A = rlnlm*(log(I*(zz-L))-Pi*I/2) + rlnlm2*(log(-I*(zz-L2))+Pi*I/2) + sfunczero;
      h = subst(if (lv, icdl[lv], ct), x, (zz-circc));
      icA[swidx] = A;
      icvals[swidx] = h;
    ,
      ctsmp++;
      ctsdeg=ctsdeg+length(if (lv, icdl[lv], icdct));
      ctsdig=ctsdig+length(if (lv, icdl[lv], icdct))*icdig;
      h = icvals[swidx] + subst(if (lv, icdl[lv], icdct), x, precision(zz-circc, icdig));
      icvals[swidx] = h;
      A = icA[swidx];
    );"""),

    # ---- per-iteration report -------------------------------------------
    ("""      print(n "=loopcnt "re" decimal digits, "ctsamples" ctsamples, "thsamples" thsamples");""",
     """      print(n "=loopcnt "re" decimal digits, "ctsamples" ctsamples, "thsamples" thsamples");
      print("CNT it=" n " re=" floor(re) " ths=" thsamples " cts=" ctsamples
            " rp=" default(realprecision) " thdig=" thdig " icdig=" icdig
            " degct=" length(ct) " degthd=" if(type(thdct)=="t_POL",length(thdct),0)
            " thfull=" thfull " icfull=" icfull
            " supf=" ctsupf " supn=" ctsupn " isf=" ctisf " isfn=" ctisfn
            " hc=" cthc " hdeg=" cthdeg " hdig=" cthdig
            " smp=" ctsmp " sdeg=" ctsdeg " sdig=" ctsdig);"""),
]


def main() -> None:
    text = SRC.read_text(encoding="utf-8", errors="replace")
    for old, new in PATCHES:
        if old not in text:
            sys.exit(f"PATCH NOT FOUND:\n{old[:200]}")
        text = text.replace(old, new, 1)
    DST.write_text(text, encoding="utf-8")
    print(f"wrote {DST}")


if __name__ == "__main__":
    main()
