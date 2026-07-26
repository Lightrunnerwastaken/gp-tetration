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
# Only the FORWARD direction Phi_{e,b} is stored; the reverse follows exactly
# from the groupoid identity Phi_{b,e}(H(theta)) = -Phi_{e,b}(theta)
# (H inverted by Newton; H' stays in [0.997, 1.003]).
# ---------------------------------------------------------------------------
import json
import os

from . import state_cache

# Where the shipped mode table lives, and -- separately -- where newly computed
# bases go.
#
# The old code did both through one path under research/reference/. That was
# wrong twice over:
#   1. research/reference/ is declared immutable (anti-gaming) in the README and
#      in METHODS, yet a plain library call wrote into it.
#   2. The path is repo-relative, so it does not exist in a non-editable
#      install: `pip install .` silently shipped a library whose advertised base
#      atlas could not find its data.
# So: the table is package data (read-only), and anything computed at runtime is
# written to the user cache next to the engine state cache.
_PKG_TABLE = os.path.join(os.path.dirname(__file__), "data", "phi_modes.json")
_REPO_TABLE = os.path.join(os.path.dirname(__file__), "..", "..",
                           "research", "reference", "phi_modes.json")

# Kept as a module attribute because callers and tests refer to it; it now names
# the READ seed, never a write target.
PHI_TABLE_PATH = _PKG_TABLE if os.path.exists(_PKG_TABLE) else _REPO_TABLE


def phi_table_write_path() -> str:
    """Writable location for on-demand mode tables (never the reference dir)."""
    return str(state_cache.cache_dir() / "phi_modes.json")

UNIVERSAL_TAIL = {
    "model": "ln A_k = lnC + gamma*ln k - 2*pi*sigma*k - sqrt(pi*c*k)",
    "sigma": 0.0789, "c": 59.78, "gamma": 9.47,
    "source": "universal mode-decay law, paper IV (papers/); pair-universal fit 2.4e-4",
}


def _read_json(path):
    if path and os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    return None


def _load_table(path=None):
    """Shipped seed table, overlaid with anything computed locally.

    An explicit `path` still wins outright -- tests and tools pass one.
    """
    if path is not None:
        loaded = _read_json(path)
        if loaded is not None:
            return loaded
    else:
        seed = _read_json(PHI_TABLE_PATH)
        local = _read_json(phi_table_write_path())
        if seed is not None or local is not None:
            table = seed or {"meta": {}, "bases": {}}
            if local:
                table.setdefault("bases", {}).update(local.get("bases", {}))
            return table
    return {"meta": {"version": 1, "dps": 60, "n_grid": 64, "k_head": 20,
                     "anchor": "e",
                     "definition": ("Phi_{e,b}(theta) = slog_e(T_b(n+theta)) - (n+theta); "
                                    "a_k = (2/N) sum_j Phi_j e^{-2pi i jk/N}, mu = grid mean"),
                     "universal_tail": UNIVERSAL_TAIL},
            "bases": {}}


def phi_modes_cached(gp, base, n_grid=64, k_head=20, path=None):
    """Complex Phi_{e,b} modes for one base, computed on demand and cached."""
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
        # Parameters that determine the entry, rather than a date: a
        # hardcoded "measured" stamp claimed every on-demand base had been
        # computed on one particular day in July 2026.
        "n_grid": n_grid,
        "k_head": k_head,
        "dps": mp.mp.dps,
    }
    # Only newly computed bases are persisted, and never into the shipped seed:
    # writing there would mutate research/reference/, which the project declares
    # immutable, and would fail outright on a read-only install.
    p = path or phi_table_write_path()
    if os.path.dirname(p):
        os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        json.dump(table, f, indent=1)
    return mu_v, modes


def phi_from_modes(mu_v, modes, theta):
    """Reconstruct Phi from mu and the mode head: mu + sum Re(a_k e^{2pi i k theta})."""
    theta = _to_mpf(theta)
    s = mp.mpf(mu_v)
    for k, a in enumerate(modes, start=1):
        s += (a * mp.exp(mp.mpc(0, 2) * mp.pi * k * theta)).real
    return s


def sexp_anchor(gp, base_b, y, n_lift=8, table_path=None, n_grid=64, k_head=20):
    """sexp_b(y) from the e-anchor and the Phi table alone (no base-b init).

    Lift the height: h = n_lift + y + Phi_{e,b}(y mod 1). The anchor supplies
    t = T_e(a) at a moderate height a, then an exact reverse peel in
    ln-coordinates: v_{j+1} = ln v_j - lnln b, which holds exactly as
    v_j = T_e(h-1-j) - lnln b above the precision threshold (mirror of the
    forward ladder; truncation stays below working precision).
    """
    y = _to_mpf(y)
    mu_v, modes = phi_modes_cached(gp, base_b, n_grid=n_grid, k_head=k_head, path=table_path)
    theta = y - mp.floor(y)
    h = n_lift + y + phi_from_modes(mu_v, modes, theta)
    ln_b = mp.log(_base_value(base_b))
    lnln_b = mp.log(ln_b)
    # Pick the switch level directly: a = h-1-j in (2.2, 3.2], then climb with
    # mpmath if needed until T_e(a) exceeds the threshold -- never an engine
    # call at a tower height that cannot be represented.
    thresh = mp.mpf(10) ** (mp.mp.dps // 2 + 10)
    j = int(mp.ceil(h - 1 - mp.mpf("3.2")))
    j = max(j, 1)
    a = h - 1 - j
    t = mp.mpf(gp.sexp("exp(1)", a).real)
    while t <= thresh and j > 1:
        t = mp.e ** t          # one e-level up: T_e(a+1)
        j -= 1
    v = t - lnln_b
    while j < n_lift:
        v = mp.log(v) - lnln_b
        j += 1
    return mp.e ** v


def _phi_deriv_from_modes(modes, theta):
    """Phi'(theta) from the mode head: sum Re(2 pi i k a_k e^{2pi i k theta})."""
    theta = _to_mpf(theta)
    s = mp.mpf(0)
    for k, a in enumerate(modes, start=1):
        s += (a * mp.mpc(0, 2) * mp.pi * k * mp.exp(mp.mpc(0, 2) * mp.pi * k * theta)).real
    return s


def slog_anchor(gp, base_b, w, table_path=None, n_grid=64, k_head=20):
    """slog_b(w) from the e-anchor and the Phi table alone (no base-b init).

    Climb the b-tower (y > 1e4, m levels), take the exact two-level entry into
    the e-world (W = lnln b + y ln b = lnln T_b(..+2)), peel in e and make ONE
    slog_e call to get s_e, the slog_e height of the tower; then solve
    H(x) = x + Phi(x) = s_e - (m+2) by Newton (H' lies in [0.997, 1.003], so
    3 steps suffice for full depth).
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
