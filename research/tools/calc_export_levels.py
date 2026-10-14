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
# exp-091: ranks 7, 8, ... -- the curves overflow float64 beyond x ~ 1.4.
# Probe at dps 16 (ranks 7-10): every base 1.5 ... 10 reaches rank 10 with
# 9-14 digits; b = 20 stops at rank 9 (no series around the fixed point of
# nonation). Bases >= 2 cost ~1 min each (Kneser even ranks, short chains),
# 1.5 ... 1.8 ~2-3 min (regular even ranks with multipliers near 1, long
# chains). Denser around eta_7 = 1.83976, where octation changes kind.
XMAX_RANK = 3.0
RANK_BASES = [1.5, 1.6, 1.7, 1.75, 1.8, 1.83, 1.85, 1.9, 2.0, 2.3, 2.718281828459045,
              3.2, 4.0, 5.0, 7.0, 10.0, 20.0, 50.0, 100.0]


def rank_step(rank: int) -> float:
    # ranks 9 and up: each value costs a step of the rank below; 1/16 halves
    # the cost, the monotone cubic in symlog space keeps ~5 digits
    return STEP if rank <= 8 else 2*STEP


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


def sample_rank(gp: FatouGP, base: str, rank: int) -> tuple[float, list]:
    """H_rank on [-(rank - 1), XMAX_RANK] at STEP: the leading values that are
    not real (even ranks leave the real axis on the left, around
    x = -(rank - 3)) are dropped, the table starts at the first real one and
    ends before the first value that is not finite."""
    x_start = -(rank - 2.0)
    step = rank_step(rank)
    n = int(round((XMAX_RANK - x_start)/step)) + 1
    # hypergridval: in order, one step of H_{rank-1} from the point one period
    # earlier instead of a full chain per point
    exprs = [f'(v->if(type(v)=="t_ERROR",{SENTINEL},v))'
             f"(iferr(alarm(20, hypergridval({rank}, {x_start}, {step}, {i + 1})), E, {SENTINEL}))"
             for i in range(n)]
    vals = []
    for v in gp.eval_batch(base, exprs):
        re, im = float(mp.re(v)), float(mp.im(v))
        ok = math.isfinite(re) and abs(re) <= 1e300 and abs(im) <= 1e-9*max(1.0, abs(re))
        vals.append(re if ok else None)
    first = next((i for i, v in enumerate(vals) if v is not None), len(vals))
    vals = vals[first:]
    if None in vals:
        vals = vals[:vals.index(None)]
    return x_start + first*step, vals


def export_ranks(args, data: dict) -> None:
    """exp-091: tables for ranks 7, 8, ... in data["levels"][str(rank)], merged
    by base. One engine per base builds the ranks in turn (each needs the one
    below); a rank that fails or carries fewer than HEX_MIN_DIGITS digits ends
    that base's tower."""
    levels = data.setdefault("levels", {})
    for b in (args.bases or RANK_BASES):
        with SmallGP(fatou_gp="fork", dps=args.dps, eval_timeout=3600.0) as gp:
            base = repr(b)
            for r in sorted(args.ranks):
                t0 = time.time()
                try:
                    d = float(gp.hyper_digits(base, r))
                    if d < HEX_MIN_DIGITS:
                        print(f"b={b} r={r}: only {d:.1f} digits -- tower ends", flush=True)
                        break
                    kind = int(mp.re(gp.eval_batch(base, [f"hykind[{r}]"])[0]))
                    t1 = time.time()
                    x0, v = sample_rank(gp, base, r)
                except Exception as exc:  # one bad base/rank must not cost the rest
                    print(f"b={b} r={r}: FEHLER {str(exc)[:200]}", flush=True)
                    break
                entry = {"b": b, "x0": x0, "step": rank_step(r), "v": v, "digits": round(d, 1), "kind": kind}
                tab = [e for e in levels.get(str(r), []) if e["b"] != b] + [entry]
                levels[str(r)] = sorted(tab, key=lambda e: e["b"])
                print(f"b={b} r={r}: {d:.1f} digits, kind {kind}, build {t1 - t0:.1f}s, "
                      f"{len(v)} pts from x = {x0} in {time.time() - t1:.1f}s", flush=True)


def write_outputs(data: dict) -> None:
    data_path = HERE / "calc_data.json"
    data_path.write_text(json.dumps(data, separators=(",", ":")), encoding="utf-8")
    template = (HERE / "calc_template.html").read_text(encoding="utf-8")
    html = template.replace("__DATA__", json.dumps(data, separators=(",", ":")))
    (HERE / "tetration_calculator.html").write_text(html, encoding="utf-8")
    print(f"wrote calc_data.json and tetration_calculator.html ({len(html)//1024} KB)")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bases", nargs="*", type=float)
    ap.add_argument("--dry-run", action="store_true", help="time only, write nothing")
    ap.add_argument("--dps", type=int, default=16)
    ap.add_argument("--merge", action="store_true",
                    help="add/replace only the given --bases in the existing tables, "
                         "with hexation for every one of them")
    ap.add_argument("--ranks", nargs="*", type=int,
                    help="exp-091: tables for these ranks (7, 8, ...) instead of pentation "
                         "and hexation; --bases defaults to RANK_BASES")
    args = ap.parse_args()
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    mp.mp.dps = 40

    data_path = HERE / "calc_data.json"
    data = json.loads(data_path.read_text(encoding="utf-8"))
    if args.ranks:
        export_ranks(args, data)
        if not args.dry_run:
            write_outputs(data)
        return
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
