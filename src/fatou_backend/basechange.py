"""Base-change ladder (M5.1): Phi_{b,d} and mu via the peel ladder.

Consumer-side implementation of the base-change phase method of the
mixed-base tetration paper series (see papers/, Papers I-II):

    Phi_{b,d}(theta) = slog_b(T_d(n+theta)) - (n+theta)

computed as: tower up y = T_d(m+theta) explicitly until y > TOWER_CUT,
exact two-level entry W = alpha_hat + c*y with c = log_b d,
alpha_hat = log_b c  (identically log_b(log_b(T_d(m+theta+2))) — no
rounding of the huge tower), numerical peeling W <- log_b W until
W <= PEEL_CUT, then one slog_b call:

    slog_b(T_d(m+theta+2)) = slog_b(W_final) + 2 + k   (k peels)
    Phi = slog_b(W_final) + k + 2 - (m + theta + 2)

mu_{b,d} is the mean of Phi over a uniform theta grid (the k=0 Fourier
mode; higher aliased modes decay fast on the regular branch).
"""
from __future__ import annotations

import mpmath as mp

TOWER_CUT = 1e4
PEEL_CUT = 50


def _to_mpf(x) -> mp.mpf:
    if isinstance(x, str):
        return mp.mpf(x)
    return mp.mpf(x)


def phi(gp_b, gp_d, base_b, base_d, theta) -> mp.mpf:
    """Phi_{b,d}(theta) via the peel ladder.

    gp_b/gp_d: FatouGP-like objects (need .sexp(base, x) and .slog(base, y)).
    base_b/base_d: GP base expressions ("exp(1)", 2, ...) matching gp_* calls.
    theta: in [0, 1).
    """
    theta = _to_mpf(theta)
    ln_b = mp.log(_base_value(base_b))
    ln_d = mp.log(_base_value(base_d))

    y = mp.mpf(gp_d.sexp(base_d, theta).real)  # T_d(theta)
    m = 0
    while y <= TOWER_CUT:
        y = mp.e ** (ln_d * y)  # d^y
        m += 1

    c = ln_d / ln_b                      # log_b d
    alpha_hat = mp.log(c) / ln_b         # log_b c
    w = alpha_hat + c * y                # == log_b(log_b(T_d(m+theta+2)))

    k = 0
    while w > PEEL_CUT:
        w = mp.log(w) / ln_b
        k += 1

    slog_w = mp.mpf(gp_b.slog(base_b, w).real)
    return slog_w + k - m - theta


def mu(gp_b, gp_d, base_b, base_d, n_grid: int = 64) -> mp.mpf:
    """mu_{b,d} = grid mean of Phi_{b,d} over theta_j = j/n_grid."""
    total = mp.mpf(0)
    for j in range(n_grid):
        total += phi(gp_b, gp_d, base_b, base_d, mp.mpf(j) / n_grid)
    return total / n_grid


def _base_value(base) -> mp.mpf:
    """Numeric value of a GP base expression."""
    if isinstance(base, str):
        s = base.strip().lower()
        if s in ("exp(1)", "e"):
            return mp.e
        return mp.mpf(base)
    return mp.mpf(base)


# ---------------------------------------------------------------------------
# M5.3: Phi mode table (JSON, on-demand + cache; k_head=20, dps=60,
# universal tail per the paper-IV mode-decay law).
# Gespeichert wird die VORWAERTS-Richtung
# Phi_{e,b}; die Rueckrichtung folgt exakt aus dem Groupoid
# Phi_{b,e}(H(theta)) = -Phi_{e,b}(theta) (H-Inversion via Newton,
# H' in [0.997, 1.003]).
# ---------------------------------------------------------------------------
import json
import os

PHI_TABLE_PATH = os.path.join(os.path.dirname(__file__), "..", "..",
                              "research", "reference", "phi_modes.json")

UNIVERSAL_TAIL = {
    "model": "ln A_k = lnC + gamma*ln k - 2*pi*sigma*k - sqrt(pi*c*k)",
    "sigma": 0.0789, "c": 59.78, "gamma": 9.47,
    "source": "universal mode-decay law, paper IV (papers/); pair-universal fit 2.4e-4",
}


def _load_table(path=None):
    path = path or PHI_TABLE_PATH
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    return {"meta": {"version": 1, "dps": 60, "n_grid": 64, "k_head": 20,
                     "anchor": "e",
                     "definition": ("Phi_{e,b}(theta) = slog_e(T_b(n+theta)) - (n+theta); "
                                    "a_k = (2/N) sum_j Phi_j e^{-2pi i jk/N}, mu = Gittermittel"),
                     "universal_tail": UNIVERSAL_TAIL},
            "bases": {}}


def phi_modes_cached(gp, base, n_grid=64, k_head=20, path=None):
    """Komplexe Phi_{e,b}-Moden fuer eine Basis, on-demand + JSON-Cache."""
    table = _load_table(path)
    key = str(base)
    if key in table["bases"]:
        e = table["bases"][key]
        mu_v = mp.mpf(e["mu"])
        modes = [mp.mpc(mp.mpf(m["re"]), mp.mpf(m["im"])) for m in e["modes"]]
        return mu_v, modes
    vals = [phi(gp, gp, "exp(1)", base, mp.mpf(j) / n_grid) for j in range(n_grid)]
    mu_v = sum(vals) / n_grid
    modes = []
    for k in range(1, k_head + 1):
        s = mp.mpc(0)
        for j, v in enumerate(vals):
            s += v * mp.exp(mp.mpc(0, -2) * mp.pi * j * k / n_grid)
        modes.append(s * 2 / n_grid)
    table["bases"][key] = {
        "mu": mp.nstr(mu_v, 45),
        "modes": [{"re": mp.nstr(m.real, 45), "im": mp.nstr(m.imag, 45)} for m in modes],
        "measured": "2026-07-13",
    }
    p = path or PHI_TABLE_PATH
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        json.dump(table, f, indent=1)
    return mu_v, modes


def phi_from_modes(mu_v, modes, theta):
    """Phi-Rekonstruktion aus mu + Moden-Kopf: mu + sum Re(a_k e^{2pi i k theta})."""
    theta = _to_mpf(theta)
    s = mp.mpf(mu_v)
    for k, a in enumerate(modes, start=1):
        s += (a * mp.exp(mp.mpc(0, 2) * mp.pi * k * theta)).real
    return s


def sexp_anchor(gp, base_b, y, n_lift=8, table_path=None, n_grid=64, k_head=20):
    """sexp_b(y) NUR aus dem e-Anker + Phi-Tabelle (kein Basis-b-Init).

    Hoehe liften: h = n_lift + y + Phi_{e,b}(y mod 1); Anker liefert
    t = T_e(a) an moderater Hoehe a; dann exakter Rueckpeel in
    ln-Koordinaten: v_{j+1} = ln v_j - lnln b, wobei oberhalb der
    Praezisionsschwelle v_j = T_e(h-1-j) - lnln b exakt gilt
    (Spiegel der Vorwaerts-Leiter; Trunkierung sub-Praezision).
    """
    y = _to_mpf(y)
    mu_v, modes = phi_modes_cached(gp, base_b, n_grid=n_grid, k_head=k_head, path=table_path)
    theta = y - mp.floor(y)
    h = n_lift + y + phi_from_modes(mu_v, modes, theta)
    ln_b = mp.log(_base_value(base_b))
    lnln_b = mp.log(ln_b)
    # Schaltlevel direkt waehlen: a = h-1-j in (2.2, 3.2], dann noetigenfalls
    # via mpmath hochklettern bis T_e(a) > Schwelle (nie Engine-Calls auf
    # unrepraesentierbaren Turmhoehen).
    thresh = mp.mpf(10) ** (mp.mp.dps // 2 + 10)
    j = int(mp.ceil(h - 1 - mp.mpf("3.2")))
    j = max(j, 1)
    a = h - 1 - j
    t = mp.mpf(gp.sexp("exp(1)", a).real)
    while t <= thresh and j > 1:
        t = mp.e ** t          # eine e-Ebene hoch: T_e(a+1)
        j -= 1
    v = t - lnln_b
    while j < n_lift:
        v = mp.log(v) - lnln_b
        j += 1
    return mp.e ** v


def _phi_deriv_from_modes(modes, theta):
    """Phi'(theta) aus dem Moden-Kopf: sum Re(2 pi i k a_k e^{2pi i k theta})."""
    theta = _to_mpf(theta)
    s = mp.mpf(0)
    for k, a in enumerate(modes, start=1):
        s += (a * mp.mpc(0, 2) * mp.pi * k * mp.exp(mp.mpc(0, 2) * mp.pi * k * theta)).real
    return s


def slog_anchor(gp, base_b, w, table_path=None, n_grid=64, k_head=20):
    """slog_b(w) NUR aus dem e-Anker + Phi-Tabelle (kein Basis-b-Init).

    In b-Tuermen hochklettern (y > 1e4, m Ebenen), exakter Zwei-Level-
    Einstieg in die e-Welt (W = lnln b + y ln b = lnln T_b(..+2)),
    e-Peel + EIN slog_e-Call -> s_e = slog_e-Hoehe des Turms; dann
    H(x) = x + Phi(x) = s_e - (m+2) per Newton loesen (H' in
    [0.997, 1.003], 3 Schritte reichen fuer Volltiefe).
    """
    w = _to_mpf(w)
    mu_v, modes = phi_modes_cached(gp, base_b, n_grid=n_grid, k_head=k_head, path=table_path)
    ln_b = mp.log(_base_value(base_b))
    y = w
    m = 0
    while y <= TOWER_CUT:
        y = mp.e ** (ln_b * y)
        m += 1
    W = mp.log(ln_b) + y * ln_b
    k = 0
    while W > PEEL_CUT:
        W = mp.log(W)
        k += 1
    s_e = mp.mpf(gp.slog("exp(1)", W).real) + 2 + k
    target = s_e - (m + 2)
    x = target - mu_v
    for _ in range(3):
        g = x + phi_from_modes(mu_v, modes, x - mp.floor(x)) - target
        x -= g / (1 + _phi_deriv_from_modes(modes, x - mp.floor(x)))
    return x
