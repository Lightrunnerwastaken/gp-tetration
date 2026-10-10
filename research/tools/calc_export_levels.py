"""Pentation and hexation curves for the calculator (exp-080).

Adds two tables to calc_data.json, one entry per base of the existing grid:
  {"b": base, "x0": first x, "step": dx, "v": [value or null, ...]}
sampled with the fork at dps 16 -- plenty for a plot whose readout shows six
digits. null marks a value beyond float64, a complex value (hexation below its
real domain) or an engine error. Then rebuilds tetration_calculator.html from
calc_template.html.

    python research/tools/calc_export_levels.py            # all grid bases
    python research/tools/calc_export_levels.py --bases 2 10 --dry-run
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import mpmath as mp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

from fatou_backend.gp_backend import FatouGP  # noqa: E402

STEP = 1/32
ETA5 = 1.7373559251516875731
# extra bases on both sides of eta5, where the hexation changes character
# (bounded below, unbounded above) and must not be interpolated across
EXTRA = [1.73, 1.737, 1.738, 1.745]
SENTINEL = "1e999"
KNMAXK = 600       # exp-081/083: centre planning makes b = 10 cost K ~ 200
HEX_BMAX = 105.0    # b = 100 works since exp-084 (series around p needs
                    # sr = 1/256; 12.1 digits at dps 16) -- the whole grid
HEX_MIN_DIGITS = 5


class SmallGP(FatouGP):
    """256 MB PARI stack instead of 2 GB, as in calc_export.py: on Windows the
    maximum is reserved up front, and 37 open workers exhausted the page file
    (first full run, 2026-10-09: "nicht genuegend Speicherressourcen" from
    b = 2.81 on)."""
    def _init_lines(self, base):
        lines = super()._init_lines(base)
        return ["default(parisizemax, 268435456);"] + [l for l in lines if "parisizemax" not in l]


def sample(gp: FatouGP, base: str, fn: str, x0: float, x1: float) -> list:
    n = int(round((x1 - x0)/STEP)) + 1
    xs = [x0 + i*STEP for i in range(n)]
    # alarm: einzelne Werte weit jenseits des Ueberlaufs rechnen sonst endlos
    # (b = 10: hex(2.28), hex(3.28) -- alle Nachbarn brechen sauber mit Fehler ab)
    # alarm() gibt bei Zeitueberschreitung ein t_ERROR zurueck statt zu werfen
    exprs = [f'(v->if(type(v)=="t_ERROR",{SENTINEL},v))'
             f"(iferr(alarm(10, {fn}({mp.nstr(mp.mpf(x), 20)})), E, {SENTINEL}))" for x in xs]
    out = []
    for v in gp.eval_batch(base, exprs):
        re, im = float(mp.re(v)), float(mp.im(v))
        if not math.isfinite(re) or abs(re) > 1e300 or abs(im) > 1e-9*max(1.0, abs(re)):
            out.append(None)
        else:
            out.append(re)
    # the curves run off to +inf: cut at the first value that is not finite.
    # Beyond it the fork once returned a zero without valid digits (b = 100:
    # pent(1.96875) = 0.E495855885102019914), which stood in the table as 0
    if None in out:
        out = out[:out.index(None)]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bases", nargs="*", type=float)
    ap.add_argument("--dry-run", action="store_true", help="time only, write nothing")
    ap.add_argument("--dps", type=int, default=16)
    ap.add_argument("--merge", action="store_true",
                    help="add/replace only the given --bases in the existing tables, "
                         "with hexation for every one of them")
    args = ap.parse_args()
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    mp.mp.dps = 40

    data_path = HERE / "calc_data.json"
    data = json.loads(data_path.read_text(encoding="utf-8"))
    bases = args.bases or sorted({e["b"] for e in data["bases"]} | set(EXTRA))

    pent, hexa = [], []
    # hexation above eta5 on every second grid base (plus the eta5 extras):
    # Kneser costs 10-100 s per base, interpolation in lnln b covers the rest
    above = [b for b in bases if ETA5 < b <= HEX_BMAX and b not in EXTRA]
    hex_bases = ({b for b in bases if b < ETA5} | set(above[::2]) | {b for b in EXTRA if b in bases}
                 | ({above[-1]} if above else set()))
    if args.merge:
        hex_bases = {b for b in bases if b <= HEX_BMAX}
    for b in bases:
        # one engine per base, closed afterwards: FatouGP keeps a worker per
        # base alive, and they piled up
        with SmallGP(fatou_gp="fork", dps=args.dps, eval_timeout=300.0) as gp:
            base = repr(b)
            t0 = time.time()
            try:
                vp = sample(gp, base, "pentx", -5.0, 4.0)
                tp = time.time() - t0
                print(f"  b={b}: pent done {tp:.1f}s", flush=True)
                pent.append({"b": b, "x0": -5.0, "step": STEP, "v": vp})
                if b not in hex_bases:
                    print(f"b={b}: pent {len(vp)} pts {tp:.1f}s, no hexation table", flush=True)
                    continue
                # Kneser nur bis KNMAXK Glieder; darueber Sheldons Naeherung, die
                # bei grossen Basen nichts taugt -- dann faellt die Basis unten raus
                gp.eval_batch(base, [f"knmaxk={KNMAXK}"])
                hd = float(gp.hexation_digits(base))
                print(f"  b={b}: hexation init done {time.time() - t0:.1f}s ({hd:.1f} digits)", flush=True)
                # ab -3: darunter wird die Hexation komplex (hex(-3) = -2 liegt
                # schon unter dem Fixpunkt der Tetration), und die ipent-Ketten
                # dort liefen im Pilot ueber 10 min fuer eine Basis
                if hd < HEX_MIN_DIGITS:
                    print(f"b={b}: hexation only {hd:.1f} digits -- dropped", flush=True)
                    continue
                vh = sample(gp, base, "hexx", -2.96875, 4.0 if b > ETA5 else 8.0)
                hexa.append({"b": b, "x0": -2.96875, "step": STEP, "v": vh, "digits": round(hd, 1)})
                print(f"b={b}: pent {len(vp)} pts {tp:.1f}s, hex {len(vh)} pts "
                      f"({hd:.1f} digits), total {time.time() - t0:.1f}s", flush=True)
            except Exception as exc:  # one bad base must not cost the whole table
                print(f"b={b}: FEHLER {exc}", flush=True)

    if args.dry_run:
        return
    if args.merge:
        def merged(old, new):
            keep = {e["b"]: e for e in old}
            keep.update({e["b"]: e for e in new})
            return [keep[b] for b in sorted(keep)]
        pent = merged(data.get("pent", []), pent)
        hexa = merged(data.get("hex", []), hexa)
    data["pent"] = pent
    data["hex"] = hexa
    data["eta5"] = ETA5
    data_path.write_text(json.dumps(data, separators=(",", ":")), encoding="utf-8")
    template = (HERE / "calc_template.html").read_text(encoding="utf-8")
    html = template.replace("__DATA__", json.dumps(data, separators=(",", ":")))
    (HERE / "tetration_calculator.html").write_text(html, encoding="utf-8")
    print(f"wrote calc_data.json and tetration_calculator.html ({len(html)//1024} KB)")


if __name__ == "__main__":
    main()
