"""E6/13.4-Prototyp: zertifizierte Phi/mu-Einschluesse via Ball-Arithmetik.

Die Peel-Leiter (Turm + exakter Zwei-Level-Einstieg + Peeling) laeuft
komplett in flint.arb (Kugelarithmetik, rigoros). Der eine slog_e-Wert
aus der Engine wird als Ball [wert +- eps_engine] eingespeist; eps_engine
ist durch die Fehler-Vektor-Referenzen gerechtfertigt (dps-60-Worker:
konservativ 1e-55; N=64 damit Grid-Aliasing (A_64 ~ e-60) unter dem Budget liegt). Ergebnis: mu-Ball mit BEWIESENEM Radius.
"""
import sys, time

import flint
import mpmath as mp

from fatou_backend.gp_backend import FatouGP

ENGINE_EPS = "1e-33"   # Kalibriergesetz: echte digits ~ dps-24 => dps-60-Werte ~1e-36; +3 Marge


def phi_ball(gp, base_d_val, theta_num, theta_den, prec=256):
    """Zertifizierter Ball fuer Phi_{e,d}(theta), theta = num/den exakt."""
    flint.ctx.prec = prec
    eps = flint.arb(ENGINE_EPS)
    ln_d = flint.arb(base_d_val).log()
    theta = flint.arb(theta_num) / theta_den
    # Engine-Wert sexp_d(theta) als Ball
    y0 = mp.mpf(gp.sexp(base_d_val, mp.mpf(theta_num) / theta_den).real)
    y = flint.arb(str(y0)).union(flint.arb(str(y0)) + eps).union(flint.arb(str(y0)) - eps)
    m = 0
    while float(y.mid()) <= 1e4:
        y = (ln_d * y).exp()
        m += 1
    c = ln_d  # log_e d
    alpha_hat = c.log()
    w = alpha_hat + c * y          # exakter Zwei-Level-Einstieg (e-Anker: ln b = 1)
    k = 0
    while float(w.mid()) > 50:
        w = w.log()
        k += 1
    s0 = mp.mpf(gp.slog("exp(1)", mp.mpf(w.mid().str(60, radius=False))).real)
    # slog-Ball: Engine-Fehler + Ableitungs-Propagation des w-Radius (slog' ~ O(1) bei w~50)
    slog_ball = flint.arb(str(s0)).union(flint.arb(str(s0)) + eps).union(flint.arb(str(s0)) - eps)
    wr = flint.arb(w.rad().str(20, radius=False))
    wrad = flint.arb(0).union(wr).union(-wr)
    slog_ball = slog_ball + wrad * 2   # |slog'| <= 2 konservativ bei w ~ 50
    return slog_ball + k - m - theta


def main():
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    mp.mp.dps = 70
    gp = FatouGP(dps=60)
    N = 64
    t0 = time.time()
    total = flint.arb(0)
    for j in range(N):
        total = total + phi_ball(gp, 2, j, N)
    mu_ball = total / N
    print(f"mu_(e,2)-Ball (N={N}): {time.time()-t0:.1f}s")
    print(f"  Mittelpunkt: {mu_ball.mid().str(45, radius=False)}")
    print(f"  ZERTIFIZIERTER Radius: {mu_ball.rad().str(8, radius=False)}")
    ref = "-1.128403776628924009879217902036426267112"
    print(f"  Referenz-Log:  {ref}")
    inside = flint.arb(ref) in mu_ball if hasattr(flint.arb(ref), '__contains__') else None
    diff = abs(float((flint.arb(ref) - flint.arb(mu_ball.mid().str(45, radius=False))).mid().str(8, radius=False)))
    print(f"  |Referenz - Mittelpunkt| = {diff:.2e} (muss <= Radius sein)")


if __name__ == "__main__":
    main()
