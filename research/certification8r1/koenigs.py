"""Certified Koenigs coordinate with a *derived* stopped-recursion remainder.

The pre-certificate stopped the backward Koenigs recursion at
``|delta_n| < 1e-16`` and charged the truncation to two asserted constants,

    gamma_error       = 128  * |delta_n| / (1 - 1/|lambda|)**3
    gamma_prime_error = 1024 * |delta_n| / (1 - 1/|lambda|)**5 ,

with no derivation of the factors ``128`` and ``1024``.  This module derives
both from quantities the run itself produces.

Notation
--------
``phi(w) = log(w)/c`` is the backward map, ``L`` its fixed point,
``phi'(L) = 1/lambda`` with ``|lambda| > 1``, ``delta_n = w_n - L``, and

    psi(d)      = phi(L+d) - L - d/lambda            (double zero at 0)
    Phi_beta(d) = d_beta [ phi_beta(L_beta+d) - L_beta ]   (simple zero at 0)

Both are evaluated as Taylor models whose constant term is inflated to the
closed disc ``|d| <= rho``, so the resulting bounds hold uniformly in ``d``
and in the base parameter at once.  Writing ``K = sup|psi|/rho**2`` and
``B = sup|Phi_beta|/rho`` the Schwarz lemma gives

    |psi(d)|  <= K |d|**2 ,     |psi'(d)| <= 2K|d| (1 + 2|d|/rho) ,
    |xi(d)|   <= B ,            |xi(d)-xi(0)| <= 4B|d|/rho ,

where ``xi(d) = Phi_beta(d)/d``.  The three telescoping sums below then bound
the stopping error of ``log kappa``, of its state derivative and of its base
derivative.  Every input to those bounds is an interval produced by the run.

Value
    delta_{m+1} = delta_m/lambda + psi(delta_m), hence
    lambda^{m+1} delta_{m+1} = lambda^m delta_m (1 + lambda eta(delta_m) delta_m)
    with |eta| <= K, so with theta = 1/|lambda| + K|delta_n| < 1,

        sum_{m>=n} |lambda eta delta_m| <= t = |lambda| K |delta_n| / (1-theta),
        |log kappa - log S^(n)| <= t/(1-t)  =: E_kappa .

State derivative
    g_m = s_m/delta_m satisfies
    log g_{m+1} - log g_m = log(1+lambda psi'(delta_m)) - log(1+lambda eta delta_m),
    so with u = |lambda| K |delta_n| (3 + 4|delta_n|/rho) / (1-theta),

        |log g - log g_n| <= u/(1-u) =: E_g ,
        |g - g_n| <= |g_n| (exp(E_g) - 1) .

Base derivative
    G_m = (v_m - dL)/delta_m and dlogS^(m) = m dloglam + G_m obey

    dlogS^(m+1) - dlogS^(m)
        = [ (xi(delta_m)-xi(0)) + dloglam eta delta_m
            + G_m (psi'(delta_m) - eta delta_m) ] / (1/lambda + eta delta_m),

    because xi(0) = -dloglam/lambda cancels the leading term.  Summing with
    |G_m| <= |G_n| + (m-n) C, C = 2|dloglam|, gives the bound implemented in
    :meth:`KoenigsBase.stopping_bounds`.

The derivative ``Gamma'`` in the state variable is *not* given a separate
telescoping analysis.  Instead the whole quotient error is computed once more
on a closed state disc of radius ``rho_y`` around the evaluation point, and
the Cauchy estimate ``|d_y f| <= sup_{disc} |f| / rho_y`` converts it.
"""

from __future__ import annotations

from flint import acb, arb

from tmodel import TM, TaylorModelError, fixed_point_exp


class KoenigsError(RuntimeError):
    """Raised when the Koenigs certificate cannot be closed."""


class KoenigsBase:
    """Base-dependent Koenigs data on one Lobatto panel, as Taylor models."""

    def __init__(self, beta: TM, order: int, psi_radius: str = "0.5"):
        self.order = order
        self.beta = beta
        self.c = beta.exp()
        self.psi_radius = arb(psi_radius)

        c_series = self.c.series()
        seed = -(-c_series).lambertw(-1) / c_series
        self.L, self.contraction = fixed_point_exp(self.c, seed, order)

        self.lam = self.c * self.L
        self.log_lam = self.lam.log()
        self.log_lam_lower = self.log_lam.lower()
        if not bool(self.log_lam_lower > 0):
            raise KoenigsError("log lambda is not separated from zero")
        self.inverse_lambda_upper = self.lam.inverse().bound()
        if not bool(self.inverse_lambda_upper < 1):
            raise KoenigsError("Koenigs multiplier is not expanding")
        self.lambda_upper = self.lam.bound()

        self.dL = self.c * self.L * self.L / (1 - self.lam)
        self.dlam = self.c * (self.dL + self.L)
        self.dlog_lam = self.dlam / self.lam
        self.dlog_lam_upper = self.dlog_lam.bound()

        self.psi_constant, self.phi_beta_constant = self._disc_constants()

    # ------------------------------------------------------- disc constants

    def _disc_disc(self) -> TM:
        """The model of ``L + D(0, rho)``."""

        return TM(list(self.L.c), self.L.r + self.psi_radius, self.order)

    def _disc_constants(self) -> tuple[arb, arb]:
        rho = self.psi_radius
        disc = self._disc_disc()
        if not bool(disc.lower() > 0):
            raise KoenigsError(
                "psi disc reaches the logarithmic branch point at zero"
            )
        offset = disc - self.L
        image = disc.log() / self.c
        psi = image - self.L - offset / self.lam
        # d_beta [ log(L+d)/c - L ] at fixed d
        #   = -log(L+d)/c + dL/(c(L+d)) - dL
        phi_beta = -image + self.dL * disc.inverse() / self.c - self.dL
        return psi.bound() / (rho * rho), phi_beta.bound() / rho

    # ------------------------------------------------------------ estimates

    def stopping_bounds(
        self, delta_upper: arb, g_upper: arb, capital_g_upper: arb
    ) -> dict:
        """Derived remainders at a stopping distance ``|delta_n| <= d``."""

        rho = self.psi_radius
        K = self.psi_constant
        B = self.phi_beta_constant
        if not bool(delta_upper < rho / 2):
            raise KoenigsError("stopping distance exceeds half the psi disc")
        theta = self.inverse_lambda_upper + K * delta_upper
        if not bool(theta < 1):
            raise KoenigsError("delta contraction is not certified")
        geometric = delta_upper / (1 - theta)

        value_series = self.lambda_upper * K * delta_upper / (1 - theta)
        if not bool(value_series < arb(1) / 2):
            raise KoenigsError("Koenigs value series is not contractive")
        kappa_error = value_series / (1 - value_series)

        shape = 3 + 4 * delta_upper / rho
        state_series = self.lambda_upper * K * delta_upper * shape / (1 - theta)
        if not bool(state_series < arb(1) / 2):
            raise KoenigsError("Koenigs state-derivative series is diverging")
        state_log_error = state_series / (1 - state_series)
        state_error = g_upper * (state_log_error.exp() - 1)

        denominator = 1 / self.lambda_upper - K * delta_upper
        if not bool(denominator > 0):
            raise KoenigsError("Koenigs base-derivative denominator vanishes")
        drift = 2 * self.dlog_lam_upper
        weighted = K * shape * delta_upper * (
            capital_g_upper / (1 - theta)
            + drift * theta / ((1 - theta) * (1 - theta))
        )
        base_error = (
            4 * B * geometric / rho
            + self.dlog_lam_upper * K * geometric * delta_upper
            + weighted
        ) / denominator

        return {
            "psi_radius": rho,
            "psi_constant": K,
            "phi_beta_constant": B,
            "theta": theta,
            "log_kappa_error": kappa_error,
            "state_log_error": state_log_error,
            "state_error": state_error,
            "base_error": base_error,
        }


class KoenigsRun:
    """One stopped backward orbit and the quantities assembled from it."""

    __slots__ = (
        "depth",
        "delta_upper",
        "phase_margin",
        "log_s",
        "dlog_s",
        "slog_derivative",
        "d_dlog_s",
        "d_slog_derivative",
        "gamma",
        "gamma_prime",
        "capital_g_upper",
        "g_upper",
        "p_lower",
    )


def _orbit(
    base: KoenigsBase,
    y: TM,
    stop: arb,
    max_depth: int,
    want_second: bool,
) -> KoenigsRun:
    """Run the stopped recursion and assemble the raw stopped quantities."""

    order = base.order
    c = base.c
    L = base.L
    dL = base.dL

    w = y
    v = TM.zero(order)
    s = TM.constant(1, order)
    u = TM.zero(order)
    s2 = TM.zero(order)

    log_delta: TM | None = None
    phase_margin: arb | None = None
    reference = base.lam.inverse().ball()

    for depth in range(1, max_depth + 1):
        inverse_cw = (c * w).inverse()
        w_next = w.log() / c
        s_next = s * inverse_cw
        v_next = v * inverse_cw - w_next
        if want_second:
            scaled = inverse_cw * inverse_cw * c
            u_next = u * inverse_cw - v * s * scaled - s_next
            s2_next = s2 * inverse_cw - s * s * scaled
        else:
            u_next = u
            s2_next = s2
        previous_delta = None if log_delta is None else w - L
        w, v, s, u, s2 = w_next, v_next, s_next, u_next, s2_next
        delta = w - L

        if log_delta is None:
            # Depth one fixes the branch: the principal logarithm, exactly as
            # in the pre-certificate.  Every later step continues it by the
            # logarithm of a ratio that stays near 1/lambda, so no winding
            # index has to be guessed.
            log_delta = delta.log()
        else:
            ratio = delta * previous_delta.inverse()
            excursion = (ratio.ball() / reference).arg().abs_upper()
            margin = arb.pi() - excursion
            if not bool(margin > 0):
                raise KoenigsError(
                    "Koenigs step ratio leaves the principal sector"
                )
            if phase_margin is None or bool(margin < phase_margin):
                phase_margin = margin
            log_delta = log_delta + ratio.log()

        delta_upper = delta.ball().abs_upper()
        if not bool(delta_upper < stop):
            continue

        inverse_delta = delta.inverse()
        run = KoenigsRun()
        run.depth = depth
        run.delta_upper = delta_upper
        run.phase_margin = phase_margin or arb.pi()
        run.log_s = base.log_lam * depth + log_delta
        capital_g = (v - dL) * inverse_delta
        run.capital_g_upper = capital_g.bound()
        run.dlog_s = base.dlog_lam * depth + capital_g
        run.slog_derivative = s * inverse_delta
        run.g_upper = run.slog_derivative.bound()
        run.p_lower = run.slog_derivative.lower()
        chi = run.log_s / base.log_lam
        numerator = run.dlog_s - chi * base.dlog_lam
        run.gamma = numerator / run.slog_derivative
        if want_second:
            run.d_dlog_s = (
                (u * delta - (v - dL) * s) * inverse_delta * inverse_delta
            )
            run.d_slog_derivative = (
                (s2 * delta - s * s) * inverse_delta * inverse_delta
            )
            d_numerator = (
                run.d_dlog_s
                - (run.slog_derivative / base.log_lam) * base.dlog_lam
            )
            run.gamma_prime = (
                d_numerator * run.slog_derivative
                - numerator * run.d_slog_derivative
            ) / (run.slog_derivative * run.slog_derivative)
        else:
            run.d_dlog_s = None
            run.d_slog_derivative = None
            run.gamma_prime = None
        return run

    raise KoenigsError("stopped Koenigs recursion did not contract")


def _gamma_error(base: KoenigsBase, run: KoenigsRun) -> tuple[arb, dict]:
    """Charge the stopped-coordinate errors to ``Gamma = N/P``."""

    certificate = base.stopping_bounds(
        run.delta_upper, run.g_upper, run.capital_g_upper
    )
    kappa_error = certificate["log_kappa_error"]
    state_error = certificate["state_error"]
    base_error = certificate["base_error"]

    chi_error = kappa_error / base.log_lam_lower
    numerator_error = base_error + base.dlog_lam_upper * chi_error
    if not bool(run.p_lower > state_error):
        raise KoenigsError("Koenigs quotient denominator is not separated")
    gamma_error = (
        numerator_error + run.gamma.bound() * state_error
    ) / (run.p_lower - state_error)
    certificate["numerator_error"] = numerator_error
    certificate["gamma_error"] = gamma_error
    return gamma_error, certificate


def koenigs_field(
    base: KoenigsBase,
    y: TM,
    stop: arb,
    max_depth: int = 500,
) -> tuple[TM, dict]:
    """Certified ``Gamma`` alone (used for the graded Hilbert trace)."""

    run = _orbit(base, y, stop, max_depth, want_second=False)
    gamma_error, certificate = _gamma_error(base, run)
    gamma = TM(run.gamma.c, run.gamma.r + gamma_error, base.order)
    trace = {
        "depth": run.depth,
        "terminal_distance": run.delta_upper,
        "phase_margin_lower": run.phase_margin,
        "inverse_lambda_upper": base.inverse_lambda_upper,
        "gamma_remainder": gamma_error,
        "gamma_model_remainder": gamma.r,
        **{
            key: certificate[key]
            for key in (
                "psi_radius",
                "psi_constant",
                "phi_beta_constant",
                "theta",
                "log_kappa_error",
                "state_error",
                "base_error",
            )
        },
    }
    return gamma, trace


def koenigs_pair(
    base: KoenigsBase,
    y: TM,
    state_radius: arb,
    stop: arb,
    max_depth: int = 500,
) -> tuple[TM, TM, dict]:
    """Certified ``(Gamma, Gamma')``.

    ``Gamma`` is taken from the thin run.  Its stopping error, and the error
    of ``Gamma'``, are taken from a second run on the closed state disc of
    radius ``state_radius`` around ``y``: the stopped-coordinate defect is
    analytic there, so the Cauchy estimate turns a uniform bound on the disc
    into a bound on the derivative.
    """

    run = _orbit(base, y, stop, max_depth, want_second=True)
    fattened = TM(list(y.c), y.r + state_radius, base.order)
    disc_run = _orbit(base, fattened, stop, max_depth, want_second=False)
    gamma_error, certificate = _gamma_error(base, run)
    disc_error, disc_certificate = _gamma_error(base, disc_run)
    gamma_prime_error = disc_error / state_radius

    gamma = TM(run.gamma.c, run.gamma.r + gamma_error, base.order)
    gamma_prime = TM(
        run.gamma_prime.c,
        run.gamma_prime.r + gamma_prime_error,
        base.order,
    )
    trace = {
        "depth": run.depth,
        "disc_depth": disc_run.depth,
        "terminal_distance": run.delta_upper,
        "disc_terminal_distance": disc_run.delta_upper,
        "phase_margin_lower": min(run.phase_margin, disc_run.phase_margin),
        "inverse_lambda_upper": base.inverse_lambda_upper,
        "state_radius": state_radius,
        "gamma_remainder": gamma_error,
        "gamma_prime_remainder": gamma_prime_error,
        "gamma_model_remainder": gamma.r,
        "gamma_prime_model_remainder": gamma_prime.r,
        "disc_gamma_error": disc_error,
        **{
            key: certificate[key]
            for key in (
                "psi_radius",
                "psi_constant",
                "phi_beta_constant",
                "theta",
                "log_kappa_error",
                "state_error",
                "base_error",
            )
        },
    }
    return gamma, gamma_prime, trace
