"""M5.3a: mu-hub raster — m(b) = mu_{e,b} over a dense base grid.

Serves both our M5 (hub table, N bases instead of N^2 pairs) and the
basechange-modes research step "m(b)-Kurve dicht rastern". Validated
against their published anchors (m(2), m(3), m(5), m(10), m(1.5)).
Output: research/reference/mu_hub.json
"""
import json, sys, time

import mpmath as mp

from fatou_backend.gp_backend import FatouGP
from fatou_backend import basechange

ANCHORS = {
    "1.5": "-9.776277436687167726082341422480654566832",
    "2":   "-1.128403776628924009879217902036426267112",
    "3":   "0.1924222937224320832082444066892813277903",
    "5":   "0.7712572961544758763804410539746626381043",
    "10":  "1.136557055713326440899375967586836969056",
}


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    mp.mp.dps = 70
    gp = FatouGP(dps=60)
    # dichtes Gitter: Kante 1.5..1.6 fein, dann log-artig bis 100
    bases = ["1.5", "1.52", "1.55", "1.58", "1.6", "1.65", "1.7", "1.8", "1.9",
             "2", "2.2", "2.5", "2.8", "3", "3.5", "4", "5", "6", "8", "10",
             "13", "16", "20", "30", "50", "100"]
    out = {}
    for b in bases:
        t0 = time.time()
        try:
            m = basechange.mu(gp, gp, "exp(1)", b, n_grid=32)
            out[b] = mp.nstr(m, 40)
            note = ""
            if b in ANCHORS:
                err = abs(m - mp.mpf(ANCHORS[b]))
                note = f"  [Anker-Abw. {mp.nstr(err, 3)}]"
            print(f"m({b}) = {mp.nstr(m, 25)}  ({time.time()-t0:.1f}s){note}", flush=True)
        except Exception as e:
            print(f"m({b}) FEHLER: {e}", flush=True)
    with open("research/reference/mu_hub.json", "w") as f:
        json.dump({"meta": {"date": "2026-07-13", "dps": 60, "n_grid": 32,
                            "method": "M5.1 peel ladder, fork 21 keeps",
                            "definition": "m(b) = mu_{e,b} (Basis d=b, Anker e)"},
                   "values": out}, f, indent=1)
    print("mu_hub.json geschrieben:", len(out), "Basen")


if __name__ == "__main__":
    main()
