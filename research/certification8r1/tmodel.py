"""Rigorous Taylor models with a certified remainder disc.

A Taylor model ``TM`` of order ``N`` represents a function ``f`` that is
analytic on a neighbourhood of the closed unit disc ``D = {|s| <= 1}`` of the
scaled panel variable ``s``.  It stores

* a coefficient list ``c[0..N-1]`` of Arb complex balls, and
* a nonnegative real radius ``r``,

with the invariant

.. math::

    f(s) \\in \\sum_{k<N} c_k s^k + D(0, r)
    \\qquad\\text{for every } |s| \\le 1 .

Every operation below returns a model whose remainder is a *proved* upper
bound, so no free residual budget is ever introduced.  The two mechanisms are

1. **exact defect bounds** for the algebraic operations.  Truncated power
   series arithmetic reproduces the true Taylor coefficients exactly, so the
   only loss is the discarded high part of a convolution, which is computed
   in full and summed; and

2. **certified Cauchy tails** for ``exp`` and ``log``.  The polynomial part
   ``P`` is entire, so for any radius ``R`` on which the outer function stays
   analytic the Taylor coefficients of ``f = g o P`` obey
   ``|a_k| <= M_R R^{-k}`` with ``M_R = sup_{|s|=R} |g(P(s))|``, whence

   .. math::

       \\sum_{k \\ge N} |a_k| \\le \\frac{M_R R^{-N}}{1 - R^{-1}} .

   ``R`` is chosen from a ladder by minimising the resulting tail, and the
   admissibility of each ``R`` is itself certified from the coefficients.

The perturbation ``r`` of the argument is always propagated by a mean-value
bound on the closed ball that contains every value of the argument.

This module is used by the ``8r1`` closure to eliminate the previously free
constant ``BETA_SERIES_REMAINDER``.
"""

from __future__ import annotations

from flint import acb, acb_poly, acb_series, arb


DEFAULT_RADIUS_LADDER = (
    arb(2),
    arb(4),
    arb(8),
    arb(16),
    arb(32),
    arb(64),
    arb(128),
    arb(512),
    arb(2048),
    arb(16384),
    arb(2) ** 20,
    arb(2) ** 30,
)


class TaylorModelError(RuntimeError):
    """Raised when a model operation cannot be certified."""


def configure(order: int, precision_bits: int) -> None:
    """Set the global Arb precision *and* the series truncation length.

    ``python-flint`` truncates every series operation at the global
    ``ctx.cap``, whose default value is ``10``.  A model built with a longer
    coefficient list is therefore silently cut back to ten coefficients
    unless ``ctx.cap`` is raised, and the discarded coefficients are then
    invisible in the serialized output.  Every entry point of the closure
    calls this function before building any model, and
    :func:`assert_full_length` re-checks the achieved length afterwards.
    """

    from flint import ctx

    ctx.prec = precision_bits
    ctx.cap = order


def assert_full_length(value: acb_series, order: int, where: str) -> acb_series:
    """Fail loudly if a series was truncated below the declared order."""

    coefficients = value.coeffs()
    if len(coefficients) < order:
        # A genuinely shorter series is only admissible when the missing
        # coefficients are exact zeros, which python-flint reports by
        # dropping them.  Distinguish the two cases via the series cap.
        if value.prec < order:
            raise TaylorModelError(
                f"{where}: series truncated to {value.prec} < {order} "
                "coefficients; call configure(order, bits) first"
            )
    return value


def _zero_acb() -> acb:
    return acb(0)


class TM:
    """Order-``N`` Taylor model on the closed unit disc with remainder disc."""

    __slots__ = ("c", "r", "n")

    def __init__(self, coefficients, remainder=None, order: int | None = None):
        coefficients = list(coefficients)
        self.n = order if order is not None else len(coefficients)
        if len(coefficients) < self.n:
            coefficients = coefficients + [
                _zero_acb() for _ in range(self.n - len(coefficients))
            ]
        self.c = [acb(value) for value in coefficients[: self.n]]
        radius = arb(0) if remainder is None else arb(remainder)
        if not radius.is_finite():
            raise TaylorModelError("non-finite Taylor-model remainder")
        if radius < 0:
            raise TaylorModelError("negative Taylor-model remainder")
        self.r = radius

    # ---------------------------------------------------------------- basics

    @staticmethod
    def constant(value, order: int) -> "TM":
        return TM([acb(value)], arb(0), order)

    @staticmethod
    def zero(order: int) -> "TM":
        return TM([], arb(0), order)

    @staticmethod
    def linear(constant_term, slope, order: int) -> "TM":
        """The model of ``constant_term + slope * s``."""

        return TM([acb(constant_term), acb(slope)], arb(0), order)

    @staticmethod
    def from_series(value: acb_series, order: int, remainder=None) -> "TM":
        return TM(value.coeffs(), remainder, order)

    def series(self) -> acb_series:
        return assert_full_length(
            acb_series(self.c, self.n), self.n, "TM.series"
        )

    def poly(self) -> acb_poly:
        return acb_poly(self.c)

    def copy(self) -> "TM":
        return TM(list(self.c), self.r, self.n)

    # ------------------------------------------------------------- envelopes

    def poly_bound(self) -> arb:
        """Upper bound of ``|P(s)|`` on ``|s| <= 1`` (coefficient 1-norm)."""

        total = arb(0)
        for coefficient in self.c:
            total += coefficient.abs_upper()
        return total

    def bound(self) -> arb:
        """Upper bound of ``|f(s)|`` on ``|s| <= 1``."""

        return self.poly_bound() + self.r

    def poly_lower(self) -> arb:
        """Lower bound of ``|P(s)|`` on ``|s| <= 1``."""

        total = self.c[0].abs_lower()
        for coefficient in self.c[1:]:
            total -= coefficient.abs_upper()
        return total

    def lower(self) -> arb:
        """Lower bound of ``|f(s)|`` on ``|s| <= 1`` (may be negative)."""

        return self.poly_lower() - self.r

    def real_upper(self) -> arb:
        """Upper bound of ``Re f(s)`` on ``|s| <= 1``."""

        total = self.c[0].real.upper()
        for coefficient in self.c[1:]:
            total += coefficient.abs_upper()
        return total + self.r

    def real_lower(self) -> arb:
        """Lower bound of ``Re f(s)`` on ``|s| <= 1``."""

        total = self.c[0].real.lower()
        for coefficient in self.c[1:]:
            total -= coefficient.abs_upper()
        return total - self.r

    def ball(self) -> acb:
        """An Arb ball containing every value ``f(s)``, ``|s| <= 1``."""

        value = self.poly()(arb("0 +/- 1"))
        return acb(
            arb(value.real.mid(), value.real.rad() + self.r),
            arb(value.imag.mid(), value.imag.rad() + self.r),
        )

    def scaled_bound(self, radius: arb) -> arb:
        """Upper bound of ``|P(s)|`` on ``|s| <= radius`` (polynomial part)."""

        total = arb(0)
        power = arb(1)
        for coefficient in self.c:
            total += coefficient.abs_upper() * power
            power = power * radius
        return total

    def scaled_tail_norm(self, radius: arb) -> arb:
        """Upper bound of ``sum_{k>=1} |c_k| radius^k``."""

        total = arb(0)
        power = radius
        for coefficient in self.c[1:]:
            total += coefficient.abs_upper() * power
            power = power * radius
        return total

    # ------------------------------------------------------------ arithmetic

    def __add__(self, other) -> "TM":
        if not isinstance(other, TM):
            return self + TM.constant(other, self.n)
        coefficients = [self.c[k] + other.c[k] for k in range(self.n)]
        return TM(coefficients, self.r + other.r, self.n)

    __radd__ = __add__

    def __neg__(self) -> "TM":
        return TM([-value for value in self.c], self.r, self.n)

    def __sub__(self, other) -> "TM":
        if not isinstance(other, TM):
            return self - TM.constant(other, self.n)
        coefficients = [self.c[k] - other.c[k] for k in range(self.n)]
        return TM(coefficients, self.r + other.r, self.n)

    def __rsub__(self, other) -> "TM":
        return TM.constant(other, self.n) - self

    def __mul__(self, other) -> "TM":
        if not isinstance(other, TM):
            scalar = acb(other)
            magnitude = scalar.abs_upper()
            return TM(
                [value * scalar for value in self.c], self.r * magnitude, self.n
            )
        product = (self.poly() * other.poly()).coeffs()
        low = product[: self.n]
        high = arb(0)
        for coefficient in product[self.n :]:
            high += coefficient.abs_upper()
        remainder = (
            high
            + self.poly_bound() * other.r
            + other.poly_bound() * self.r
            + self.r * other.r
        )
        return TM(low, remainder, self.n)

    __rmul__ = __mul__

    def inverse(self) -> "TM":
        """Certified model of ``1 / f``."""

        lower = self.poly_lower()
        if not bool(lower > 0):
            raise TaylorModelError(
                "polynomial part is not bounded away from zero on the disc"
            )
        if not bool(lower - self.r > 0):
            raise TaylorModelError(
                "remainder disc reaches zero; inverse is not certified"
            )
        approximate = assert_full_length(
            acb_series([acb(1)], self.n) / self.series(), self.n, "TM.inverse"
        )
        candidate = TM.from_series(approximate, self.n)
        product = (self.poly() * candidate.poly()).coeffs()
        defect = (product[0] - acb(1)).abs_upper()
        for coefficient in product[1:]:
            defect += coefficient.abs_upper()
        if not bool(defect < 1):
            raise TaylorModelError("series inverse defect is not contractive")
        polynomial_error = candidate.poly_bound() * defect / (1 - defect)
        argument_error = self.r / (lower * (lower - self.r))
        return TM(candidate.c, polynomial_error + argument_error, self.n)

    def __truediv__(self, other) -> "TM":
        if not isinstance(other, TM):
            return self * (acb(1) / acb(other))
        return self * other.inverse()

    def __rtruediv__(self, other) -> "TM":
        return TM.constant(other, self.n) * self.inverse()

    def __pow__(self, exponent: int) -> "TM":
        if exponent < 0:
            return (self**-exponent).inverse()
        result = TM.constant(1, self.n)
        base = self
        while exponent:
            if exponent & 1:
                result = result * base
            exponent >>= 1
            if exponent:
                base = base * base
        return result

    # -------------------------------------------------- elementary functions

    def _cauchy_tail(self, modulus_at_radius, ladder=DEFAULT_RADIUS_LADDER) -> arb:
        """Minimal certified Cauchy tail ``sum_{k>=N} |a_k|``.

        ``modulus_at_radius(R)`` must return either ``None`` (radius not
        admissible) or a proved upper bound of ``sup_{|s|=R} |g(P(s))|``.
        """

        best: arb | None = None
        for radius in ladder:
            modulus = modulus_at_radius(radius)
            if modulus is None:
                continue
            tail = modulus / (radius**self.n) / (1 - 1 / radius)
            if not tail.is_finite():
                continue
            if best is None or bool(tail < best):
                best = tail
        if best is None:
            raise TaylorModelError("no admissible Cauchy radius for the tail")
        return best

    def exp(self) -> "TM":
        """Certified model of ``exp(f)``."""

        def modulus_at_radius(radius: arb):
            growth = self.scaled_tail_norm(radius)
            if not growth.is_finite():
                return None
            exponent = self.c[0].real.upper() + growth
            if not exponent.is_finite():
                return None
            modulus = exponent.exp()
            # The tower of the complex tail reaches exponents of order 1e7,
            # so no ceiling may be imposed here: Arb represents the magnitude
            # exactly, and the bound stays *relatively* tight because the
            # enclosed value is equally large.
            return modulus if modulus.is_finite() else None

        tail = self._cauchy_tail(modulus_at_radius)
        approximate = assert_full_length(
            self.series().exp(), self.n, "TM.exp"
        )
        lipschitz = (self.real_upper()).exp()
        return TM(
            approximate.coeffs(), tail + self.r * lipschitz, self.n
        )

    def log(self) -> "TM":
        """Certified model of ``log(f)`` on the branch through ``c_0``.

        The constant term is first rotated onto the positive real axis.  This
        is mathematically the identity ``log f = log(f e^{-i t}) + i t`` for a
        real ``t``, but it moves the enclosure away from the principal cut, so
        a ball that straddles the negative reals -- which happens on every
        backward Koenigs orbit, since the orbit lands on the negative axis
        before it spirals to the fixed point -- still gets a tight and, more
        importantly, *branch-consistent* enclosure.  Choosing the rotation
        from the midpoint of ``c_0`` selects the analytic continuation through
        the centre of the model, which is the continuation of the real-channel
        branch when the model is centred on the real channel.
        """

        centre = acb(
            arb(self.c[0].real.mid()), arb(self.c[0].imag.mid())
        )
        if not bool(centre.abs_lower() > 0):
            raise TaylorModelError("log argument midpoint is zero")
        angle = centre.arg()
        cosine = (-angle).cos()
        sine = (-angle).sin()
        rotation = acb(cosine, sine)
        # The exact rotation is the unit vector at -angle; the computed one
        # differs from it by at most ``offset``, and log of the exact vector
        # is exactly -i*angle for |angle| <= pi.  Charge the difference.
        offset = (
            arb(cosine.rad()) ** 2 + arb(sine.rad()) ** 2
        ).sqrt()
        modulus_lower = rotation.abs_lower()
        if not bool(modulus_lower > offset):
            raise TaylorModelError("log rotation is not separated from zero")
        rotation_defect = offset / (modulus_lower - offset)
        rotated = self * rotation
        result = rotated._principal_log() + acb(0, angle)
        return TM(result.c, result.r + rotation_defect, self.n)

    def _principal_log(self) -> "TM":
        lower = self.poly_lower()
        if not bool(lower > 0):
            raise TaylorModelError("log argument is not zero-free on the disc")
        if not bool(lower - self.r > 0):
            raise TaylorModelError("log remainder disc reaches zero")
        constant_modulus_lower = self.c[0].abs_lower()
        log_constant = self.c[0].abs_upper().log().abs_upper()
        log_constant = max(
            log_constant, constant_modulus_lower.log().abs_upper()
        )
        base_modulus = log_constant + arb.pi()

        def modulus_at_radius(radius: arb):
            growth = self.scaled_tail_norm(radius)
            if not growth.is_finite():
                return None
            ratio = growth / constant_modulus_lower
            if not bool(ratio < 1):
                return None
            return base_modulus + (-(1 - ratio).log())

        tail = self._cauchy_tail(modulus_at_radius)
        approximate = assert_full_length(
            self.series().log(), self.n, "TM.log"
        )
        argument_error = self.r / (lower - self.r)
        return TM(approximate.coeffs(), tail + argument_error, self.n)

    # ----------------------------------------------------------- diagnostics

    def certified(self) -> bool:
        return bool(self.r.is_finite()) and all(
            bool(value.is_finite()) for value in self.c
        )


def fixed_point_exp(c: TM, seed: acb_series, order: int) -> tuple[TM, arb]:
    """Certified model of the repelling fixed point ``L = exp(c L)``.

    ``seed`` supplies the polynomial part (obtained from the Lambert ``W``
    branch).  The certificate uses the *inverse* map ``phi(z) = log(z) / c``,
    which is a contraction near the repelling fixed point with factor
    ``1/|lambda|``.  If ``phi`` moves the seed by at most ``d`` and is a
    ``kappa``-contraction on the closed disc of radius ``d/(1-kappa)`` around
    it, that disc contains the unique fixed point.

    Returns the model and the certified contraction factor.
    """

    candidate = TM.from_series(seed, order)
    image = candidate.log() / c
    displacement = (image - candidate).bound()
    product = c * candidate
    contraction_lower = product.lower()
    if not bool(contraction_lower > 1):
        raise TaylorModelError("fixed-point map is not a contraction")
    kappa = 1 / contraction_lower
    enclosure = displacement / (1 - kappa)
    inflated = TM(candidate.c, candidate.r + enclosure, order)
    # Re-verify the contraction on the inflated disc.
    product_inflated = c * inflated
    if not bool(product_inflated.lower() > 1):
        raise TaylorModelError("contraction fails on the enclosing disc")
    return inflated, kappa
