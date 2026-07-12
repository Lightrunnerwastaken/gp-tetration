"""Base-change ladder (M5.1): Phi_{b,d} and mu via the peel ladder.

Consumer-side implementation of the basechange-modes research method
(Tetration/Research/basechange-modes, RESEARCH_LOG.md "Methode"):

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
