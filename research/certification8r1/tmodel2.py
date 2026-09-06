"""Two-variable Taylor models: a model in ``tau`` whose coefficients are
one-variable models in ``s``.

Motivation
----------
The base dependence must be carried as a Taylor model: a plain interval for
``beta`` is destroyed by the wrapping effect within twenty steps of the
backward Koenigs orbit, however small the interval is.  Bounding the Koenigs
field on a *complex region* -- which is what a certified Gauss--Legendre
remainder needs -- requires a second modelled variable for exactly the same
reason.  This module supplies that second variable.

A :class:`TM2` represents

.. math::

    f(s,\\tau) \\in \\sum_{j<n_\\tau} A_j(s)\\,\\tau^j + D(0,r),
    \\qquad |s| \\le 1,\\ |\\tau| \\le 1,

where each ``A_j`` is a certified one-variable :class:`~tmodel.TM`.  All the
scalar work -- inverse, ``exp``, ``log`` of the ``tau``-constant term -- is
delegated to the already-tested one-variable code, so the new material here
is only the ``tau`` recursion and the ``tau`` tail bounds.

``tau``-recursions used below (``A`` given, ``B`` sought):

* inverse   ``B_0 = A_0^{-1}``, ``B_j = -B_0 sum_{k=1..j} A_k B_{j-k}``;
* exp       ``E_0 = exp(A_0)``, ``j E_j = sum_{k=1..j} k A_k E_{j-k}``;
* log       ``L_0 = log(A_0)``,
            ``j L_j = j A_j/A_0 - (1/A_0) sum_{k=1..j-1} k L_k A_{j-k}``.

Each is the standard truncated-series identity, so the computed coefficients
are the exact Taylor coefficients of the result and the only loss is the
discarded ``tau``-tail, which is bounded by a Cauchy estimate on a certified
radius exactly as in the one-variable case.
"""

from __future__ import annotations

from flint import acb, arb

from tmodel import TM, TaylorModelError


TAU_RADIUS_LADDER = (
    arb(2),
    arb(3),
    arb(4),
    arb(6),
    arb(8),
    arb(12),
    arb(16),
    arb(32),
    arb(64),
    arb(256),
)


class TM2:
    """Order ``(n_s, n_tau)`` Taylor model on the closed unit bidisc."""

    __slots__ = ("a", "r", "ns", "nt")

    def __init__(self, coefficients: list[TM], remainder=None, order_tau=None):
        self.nt = order_tau if order_tau is not None else len(coefficients)
        if not coefficients:
            raise TaylorModelError("TM2 needs at least one coefficient")
        self.ns = coefficients[0].n
        padded = list(coefficients[: self.nt])
        while len(padded) < self.nt:
            padded.append(TM.zero(self.ns))
        self.a = padded
        radius = arb(0) if remainder is None else arb(remainder)
        if not radius.is_finite() or radius < 0:
            raise TaylorModelError("invalid TM2 remainder")
        self.r = radius

    # ---------------------------------------------------------------- basics

    @staticmethod
    def constant(value: TM, order_tau: int) -> "TM2":
        return TM2([value], arb(0), order_tau)

    @staticmethod
    def zero(order_s: int, order_tau: int) -> "TM2":
        return TM2([TM.zero(order_s)], arb(0), order_tau)

    @staticmethod
    def linear(constant_term: TM, slope: TM, order_tau: int) -> "TM2":
        return TM2([constant_term, slope], arb(0), order_tau)

    def copy(self) -> "TM2":
        return TM2([value.copy() for value in self.a], self.r, self.nt)

    # ------------------------------------------------------------- envelopes

    def poly_bound(self) -> arb:
        total = arb(0)
        for value in self.a:
            total += value.bound()
        return total

    def bound(self) -> arb:
        return self.poly_bound() + self.r

    def poly_lower(self) -> arb:
        total = self.a[0].lower()
        for value in self.a[1:]:
            total -= value.bound()
        return total

    def lower(self) -> arb:
        return self.poly_lower() - self.r

    def scaled_tail_norm(self, radius: arb) -> arb:
        total = arb(0)
        power = radius
        for value in self.a[1:]:
            total += value.bound() * power
            power = power * radius
        return total

    def ball(self) -> acb:
        """An Arb ball containing every value on the closed bidisc."""

        head = self.a[0].ball()
        spread = self.r
        for value in self.a[1:]:
            spread += value.bound()
        return acb(
            arb(head.real.mid(), head.real.rad() + spread),
            arb(head.imag.mid(), head.imag.rad() + spread),
        )

    def real_upper(self) -> arb:
        total = self.a[0].real_upper()
        for value in self.a[1:]:
            total += value.bound()
        return total + self.r

    # ------------------------------------------------------------ arithmetic

    def __add__(self, other) -> "TM2":
        if isinstance(other, TM2):
            return TM2(
                [self.a[j] + other.a[j] for j in range(self.nt)],
                self.r + other.r,
                self.nt,
            )
        if isinstance(other, TM):
            coefficients = list(self.a)
            coefficients[0] = coefficients[0] + other
            return TM2(coefficients, self.r, self.nt)
        return self + TM.constant(other, self.ns)

    __radd__ = __add__

    def __neg__(self) -> "TM2":
        return TM2([-value for value in self.a], self.r, self.nt)

    def __sub__(self, other) -> "TM2":
        return self + (-_lift(other, self))

    def __rsub__(self, other) -> "TM2":
        return (-self) + other

    def __mul__(self, other) -> "TM2":
        if not isinstance(other, TM2):
            if isinstance(other, TM):
                return TM2(
                    [value * other for value in self.a],
                    self.r * other.bound(),
                    self.nt,
                )
            scalar = acb(other)
            return TM2(
                [value * scalar for value in self.a],
                self.r * scalar.abs_upper(),
                self.nt,
            )
        low = [TM.zero(self.ns) for _ in range(self.nt)]
        tail = arb(0)
        for i, left in enumerate(self.a):
            for j, right in enumerate(other.a):
                product = left * right
                if i + j < self.nt:
                    low[i + j] = low[i + j] + product
                else:
                    tail += product.bound()
        remainder = (
            tail
            + self.poly_bound() * other.r
            + other.poly_bound() * self.r
            + self.r * other.r
        )
        return TM2(low, remainder, self.nt)

    __rmul__ = __mul__

    # ------------------------------------------------------ tau tail bounds

    def _tau_tail(self, modulus_at_radius) -> arb:
        best: arb | None = None
        for radius in TAU_RADIUS_LADDER:
            modulus = modulus_at_radius(radius)
            if modulus is None:
                continue
            tail = modulus / (radius**self.nt) / (1 - 1 / radius)
            if not tail.is_finite():
                continue
            if best is None or bool(tail < best):
                best = tail
        if best is None:
            raise TaylorModelError("no admissible tau radius for the tail")
        return best

    # ------------------------------------------------------------- inverse

    def inverse(self) -> "TM2":
        lower = self.poly_lower()
        if not bool(lower > 0):
            raise TaylorModelError("TM2 polynomial part reaches zero")
        if not bool(lower - self.r > 0):
            raise TaylorModelError("TM2 remainder reaches zero")
        head = self.a[0].inverse()
        coefficients = [head]
        for j in range(1, self.nt):
            total = TM.zero(self.ns)
            for k in range(1, j + 1):
                total = total + self.a[k] * coefficients[j - k]
            coefficients.append(-(head * total))
        candidate = TM2(coefficients, arb(0), self.nt)

        # A-posteriori defect ``|1 - A B|`` on the closed bidisc: the low part
        # is the rounding and one-variable truncation residue, the high part
        # is the discarded tau-tail of the product.  Both are computed, not
        # estimated.
        low = [TM.zero(self.ns) for _ in range(self.nt)]
        total_defect = arb(0)
        for i, left in enumerate(self.a):
            for j, right in enumerate(candidate.a):
                term = left * right
                if i + j < self.nt:
                    low[i + j] = low[i + j] + term
                else:
                    total_defect += term.bound()
        total_defect += (low[0] - 1).bound()
        for value in low[1:]:
            total_defect += value.bound()
        if not bool(total_defect < 1):
            raise TaylorModelError("TM2 inverse defect is not contractive")
        polynomial_error = (
            candidate.poly_bound() * total_defect / (1 - total_defect)
        )
        argument_error = self.r / (lower * (lower - self.r))
        return TM2(
            candidate.a, polynomial_error + argument_error, self.nt
        )

    def __truediv__(self, other) -> "TM2":
        if isinstance(other, TM2):
            return self * other.inverse()
        if isinstance(other, TM):
            return self * other.inverse()
        return self * (acb(1) / acb(other))

    def __rtruediv__(self, other) -> "TM2":
        return _lift(other, self) * self.inverse()

    # -------------------------------------------------- elementary functions

    def exp(self) -> "TM2":
        def modulus_at_radius(radius: arb):
            growth = self.scaled_tail_norm(radius)
            if not growth.is_finite():
                return None
            exponent = self.a[0].real_upper() + growth
            if not exponent.is_finite() or bool(exponent > 4096):
                return None
            return exponent.exp()

        tail = self._tau_tail(modulus_at_radius)
        coefficients = [self.a[0].exp()]
        for j in range(1, self.nt):
            total = TM.zero(self.ns)
            for k in range(1, j + 1):
                total = total + self.a[k] * coefficients[j - k] * k
            coefficients.append(total * (arb(1) / j))
        lipschitz = self.real_upper().exp()
        return TM2(coefficients, tail + self.r * lipschitz, self.nt)

    def log(self) -> "TM2":
        lower = self.poly_lower()
        if not bool(lower > 0):
            raise TaylorModelError("TM2 log argument reaches zero")
        if not bool(lower - self.r > 0):
            raise TaylorModelError("TM2 log remainder reaches zero")
        head_lower = self.a[0].lower()
        head_upper = self.a[0].bound()
        base_modulus = max(
            head_upper.log().abs_upper(), head_lower.log().abs_upper()
        ) + arb.pi()

        def modulus_at_radius(radius: arb):
            growth = self.scaled_tail_norm(radius)
            if not growth.is_finite():
                return None
            ratio = growth / head_lower
            if not bool(ratio < 1):
                return None
            return base_modulus + (-(1 - ratio).log())

        tail = self._tau_tail(modulus_at_radius)
        head_inverse = self.a[0].inverse()
        coefficients = [self.a[0].log()]
        for j in range(1, self.nt):
            total = self.a[j] * j
            for k in range(1, j):
                total = total - coefficients[k] * self.a[j - k] * k
            coefficients.append(head_inverse * total * (arb(1) / j))
        argument_error = self.r / (lower - self.r)
        return TM2(coefficients, tail + argument_error, self.nt)


def sum_tm(values: list[TM], order_s: int) -> TM:
    total = TM.zero(order_s)
    for value in values:
        total = total + value
    return total


def _lift(value, model: "TM2") -> "TM2":
    if isinstance(value, TM2):
        return value
    if isinstance(value, TM):
        return TM2.constant(value, model.nt)
    return TM2.constant(TM.constant(value, model.ns), model.nt)
