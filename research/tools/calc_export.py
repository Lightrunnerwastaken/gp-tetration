"""Export fuer den Tetration-Graphing-Calculator (statisches HTML).

Sammelt: (1) Anker-Kurve sexp_e auf dichtem Gitter, (2) mu + Phi-Moden
(k<=10) fuer ein Basis-Gitter [1.5, 100]. Schreibt calc_data.json.
Kleiner parisizemax (256MB) gegen die Pagefile-Knappheit.
"""
import json, sys, time
from pathlib import Path

import mpmath as mp

from fatou_backend.gp_backend import FatouGP
from fatou_backend import basechange


class SmallGP(FatouGP):
    def _init_lines(self, base):
        lines = super()._init_lines(base)
        return ["default(parisizemax, 268435456);"] + [l for l in lines if "parisizemax" not in l]


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    mp.mp.dps = 40
    gp = SmallGP(dps=60)
    # 1) Anker-Kurve sexp_e: y in [-2, 3.5], Schritt 1/512
    t0 = time.time()
    y0, y1, step = -2.0, 3.5, 1.0/512
    n = int(round((y1 - y0)/step)) + 1
    anchor = []
    for i in range(n):
        y = y0 + i*step
        anchor.append(float(mp.mpf(gp.sexp("exp(1)", mp.mpf(y)).real)))
    print(f"Anker-Kurve: {n} Punkte, {time.time()-t0:.0f}s")
    # 2) Basis-Gitter: Kante dicht, log bis 100
    import math
    edge = [round(1.5 + 0.02*i, 2) for i in range(15)]          # 1.50..1.78
    mid  = [round(1.8 + 0.05*i, 2) for i in range(14)]          # 1.80..2.45
    # geometrisch in lnln b von 2.5 bis 100
    s0, s1 = math.log(math.log(2.5)), math.log(math.log(100.0))
    geo = [round(math.exp(math.exp(s0 + (s1-s0)*i/40)), 4) for i in range(41)]
    bases = edge + mid + geo
    table = []
    for b in bases:
        t0 = time.time()
        try:
            vals = [basechange.phi(gp, gp, "exp(1)", repr(b), mp.mpf(j)/64) for j in range(64)]
            mu = sum(vals)/64
            modes = []
            for k in range(1, 17):
                s = mp.mpc(0)
                for j, v in enumerate(vals):
                    s += v * mp.exp(mp.mpc(0, -2)*mp.pi*j*k/64)
                s *= mp.mpf(2)/64
                modes.append([float(s.real), float(s.imag)])
            table.append({"b": b, "mu": float(mu), "modes": modes})
            print(f"  b={b}: ok ({time.time()-t0:.1f}s)", flush=True)
        except Exception as e:
            print(f"  b={b}: FEHLER {e}", flush=True)
    out = {"anchor": {"y0": y0, "step": step, "values": anchor},
           "bases": table}
    Path("research/tools/calc_data.json").write_text(json.dumps(out))
    print("calc_data.json geschrieben:", len(table), "Basen")


if __name__ == "__main__":
    main()
