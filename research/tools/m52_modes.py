"""M5.2 minimal: Phi_{b,d} mode extraction via the peel ladder (basechange.py).

Usage: python research/tools/m52_modes.py  — (e,2) demo grid, prints mu + A_k.
Validated 2026-07-13: mu matches research log to 24 digits @dps60/N64 (0.2s
warm); A2/A1, A3/A1 within ~2% of the published diagonal-limit ratios
(expected pair dependence).
"""
import sys, time

import mpmath as mp

from fatou_backend.gp_backend import FatouGP
from fatou_backend import basechange


def phi_modes(gp, base_b, base_d, n_grid=64, k_max=5):
    vals = [basechange.phi(gp, gp, base_b, base_d, mp.mpf(j)/n_grid)
            for j in range(n_grid)]
    mu = sum(vals)/n_grid
    modes = []
    for k in range(1, k_max + 1):
        s = mp.mpc(0)
        for j, v in enumerate(vals):
            s += v * mp.exp(mp.mpc(0, -2)*mp.pi*j*k/n_grid)
        modes.append(abs(s)*2/n_grid)
    return mu, modes


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    mp.mp.dps = 70
    gp = FatouGP(dps=60)
    t0 = time.time()
    mu, A = phi_modes(gp, "exp(1)", 2)
    print(f"(e,2) N=64 @dps60: {time.time()-t0:.1f}s")
    print("mu =", mp.nstr(mu, 25))
    print("A  =", [mp.nstr(a, 8) for a in A])
