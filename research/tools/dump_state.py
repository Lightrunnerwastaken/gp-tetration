"""F0: freeze a real sexpinit state so evaluation primitives can be benchmarked
on the actual data instead of on synthetic point clouds.

E1b's decisive failure mode was the point DISTRIBUTION (|w| spread over the
whole disk), which no synthetic benchmark reproduces. So the dump contains:

  ct        the Taylor correction polynomial (the thing being evaluated)
  swz/swn   staylor's cached walk endpoints and step counts -- the actual
            evaluation points, plus the arc structure (distinct step counts)
  thsw      thtaylor's cached superf points -- the second evaluation batch
  circc, circr, r, icdig, thdig, precis, samples

Numbers are written as exact integers trunc(x * 10^D) so nothing is lost or
misparsed; the loader divides by 10^D in mpmath.

    python research/tools/dump_state.py --dps 300 --out research/tools/state_300.txt
"""
from __future__ import annotations

import argparse
import subprocess
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
FORK = REPO / "src" / "fatou_backend" / "vendor" / "fatou_fork.gp"
GP = r"C:\Program Files (x86)\Pari64-2-17-3\gp.exe"

SCRIPT = r"""
default(parisizemax, 8589934592);
default(realprecision, %(dps)d);
read("%(fork)s");
default(realprecision, %(dps)d);
quietmode=1;
dbase = %(base)s;
gettime();
sexpinit(dbase,0,0,0);
dms = gettime();
DD = %(dumpdig)d;
fn = "%(out)s";
write(fn, "# fatou state dump");
write(fn, "D ", DD);
write(fn, "dps ", precis);
write(fn, "init_ms ", dms);
write(fn, "icdig ", icdig);
write(fn, "thdig ", thdig);
write(fn, "ctr ", truncate(ctr*10^30));
ei(z) = Str(truncate(real(z)*10^DD), " ", truncate(imag(z)*10^DD));
write(fn, "circc ", ei(circc));
write(fn, "circr ", ei(circr));
write(fn, "sampr ", ei(circr*ctr));
cv = Vecrev(ct);
write(fn, "ct_len ", #cv);
for (i=1, #cv, write(fn, "ct ", ei(cv[i])));
write(fn, "sw_len ", #swz);
for (i=1, #swz, write(fn, "sw ", swn[i], " ", swvalid[i], " ", ei(swz[i])));
if (type(thsw)=="t_VEC",
  write(fn, "th_len ", #thsw);
  for (i=1, #thsw, write(fn, "th ", ei(thsw[i]))),
  write(fn, "th_len 0"));
write(fn, "anchor ", ei(sexp(0.5)));
quit;
"""


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dps", type=int, default=300)
    ap.add_argument("--base", default="exp(1)")
    ap.add_argument("--dumpdig", type=int, default=0,
                    help="digits kept per number (default dps+40)")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    out = Path(args.out) if args.out else REPO / "research" / "tools" / f"state_{args.dps}.txt"
    out.unlink(missing_ok=True)
    dumpdig = args.dumpdig or args.dps + 40
    script = SCRIPT % {"dps": args.dps, "fork": FORK.as_posix(), "base": args.base,
                       "dumpdig": dumpdig, "out": out.as_posix()}
    t0 = time.time()
    proc = subprocess.run([GP, "-q", "-f"], input=script, capture_output=True,
                          text=True, encoding="utf-8", errors="replace", cwd=str(REPO))
    print(f"wall {time.time()-t0:.1f}s  rc={proc.returncode}")
    if proc.stdout.strip():
        print(proc.stdout[-800:])
    if proc.stderr.strip():
        print("STDERR", proc.stderr[-800:])
    if out.exists():
        print(f"wrote {out}  ({out.stat().st_size/1e6:.1f} MB)")
    else:
        raise SystemExit("no dump written")


if __name__ == "__main__":
    main()
