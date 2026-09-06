"""Certified complex tail ``S`` with a *proved* truncation rule.

The pre-certificate truncated the reciprocal-product tail

.. math::

    S = \\frac{y}{y'} + \\sum_{j\\ge1} \\frac1{c\\,y_j'},
    \\qquad y_{j+1}=e^{cy_j},\\quad y_{j+1}'=c\\,y_{j+1}\\,y_j',

and charged the remainder to the asserted rules ``R <= 2|a_N|`` for the value
and ``R <= 4 * 2|a_N|`` for the two partial derivatives.  Neither factor was
derived.  This module proves them.

Ratio identity
--------------
With ``a_j = 1/(c y_j')`` one has *exactly*

.. math::

    \\frac{a_{j+1}}{a_j} = \\frac{y_j'}{y_{j+1}'} = \\frac1{c\\,y_{j+1}} .

So the tail is dominated by a geometric series as soon as ``|c y_j| >= 2``
for every ``j > N``, and then

.. math::

    \\Bigl|\\sum_{j>N} a_j\\Bigr|
    \\le |a_N| \\sum_{k\\ge1} 2^{-k} = |a_N| ,

which is half of the rule that was assumed.

Forward invariance
------------------
The hypothesis is certified once, at the truncation index, by the following
elementary induction on the real channel.  The states of the segment are real
for real base and real ``x``; the Taylor models therefore enclose *real*
quantities, and a real enclosure suffices.

    Let ``c > 0`` and ``y_N`` real with ``c y_N >= log(2/c)``.  Then
    ``y_{N+1} = exp(c y_N) >= 2/c`` is real, and ``c y_{N+1} >= 2 >= log(2/c)``
    whenever ``c <= 2``, so the hypothesis reproduces itself.  Consequently
    ``y_j`` is real, increasing and ``>= 2/c`` for every ``j > N``.

Both premises -- ``c y_N >= log(2/c)`` and ``c <= 2`` -- are checked as
interval statements by :func:`certified_tail`, which refuses to return a
truncation otherwise.  On this segment ``c = exp(beta) < 0.73`` and the
tower has already passed ``10^{40000}`` at the truncation index, so the
premises hold with an enormous margin; the point is that they are *checked*.

Derivatives
-----------
Differentiating ``a_j`` with respect to the two seeds gives

.. math::

    \\partial_y a_j = -\\frac{\\partial_y y_j'}{c (y_j')^2},
    \\qquad
    \\partial_{y'} a_j = -\\frac{\\partial_{y'} y_j'}{c (y_j')^2} .

For the ``y'`` seed the tower does not depend on ``y'`` at all
(``\\partial_{y'} y_j = 0``), so ``\\partial_{y'} y'_{j+1} = c y_{j+1}
\\partial_{y'} y'_j`` and the sensitivity ratio is *exactly* the value ratio
``1/(c y_{j+1})``.

For the ``y`` seed it is **not**.  With ``\\partial_y y_{j+1} = c y_{j+1}
\\partial_y y_j`` and ``\\partial_y y'_{j+1} = c y_{j+1} [ c y'_j \\partial_y
y_j + \\partial_y y'_j ]`` one finds

.. math::

    \\frac{\\partial_y a_{j+1}}{\\partial_y a_j}
    = \\frac{1}{c\\,y_{j+1}}
      \\Bigl( 1 + \\frac{c\\,y'_j\\,\\partial_y y_j}{\\partial_y y'_j} \\Bigr)
    = \\frac{1 + c\\,y_j\\,\\theta_{j-1}}{c\\,y_{j+1}},
    \\qquad
    \\theta_{j-1} = \\frac{c y'_{j-1} \\partial_y y_{j-1}}
                        {c y'_{j-1} \\partial_y y_{j-1} + \\partial_y y'_{j-1}}
    \\in [0, 1] ,

on the real channel, where every factor is positive.  The extra factor is
``1 + c y_j theta``, up to ``1 + c y_j`` -- about twenty at the truncation
index of this segment.  An earlier revision of this module claimed the
derivative tails were dominated by the *same* geometric series as the value;
that is false for ``\\partial_y``, and it was caught by comparing the closure's
enclosures with the constant producer's, which sums one more tower level
explicitly (``test_tails.py`` keeps that comparison as a regression guard).

With ``u_j = c y_j`` the ratio is ``g(u_j) = (1 + u_j)/(c e^{u_j})`` at worst,
and ``(1+u) e^{-u}`` decreases for ``u > 0``.  The tower is monotone from the
truncation index on -- ``c e^u - u`` is increasing for ``u >= log(2/c)`` and
equals ``2 - log(2/c)`` there, which is positive when ``log(2/c) < 2`` -- so
``u_j >= u_N`` for every ``j >= N`` and

.. math::

    \\Bigl| \\sum_{j>N} \\partial_y a_j \\Bigr|
    \\le |\\partial_y a_N| \\, \\frac{r_y}{1 - r_y},
    \\qquad
    r_y = \\frac{1 + \\underline{u_N}}{c\\, e^{\\underline{u_N}}} ,

with ``\\underline{u_N}`` the certified lower bound of ``c y_N`` on the panel.
The premise ``log(2/c) < 2`` is checked next to the two invariance premises.
"""

from __future__ import annotations

from flint import acb, arb

from tmodel import TM, TaylorModelError


class TailError(RuntimeError):
    """Raised when the tail truncation cannot be certified."""


def certified_tail(
    c: TM,
    y: TM,
    y_prime: TM,
    stop: arb,
    max_terms: int = 64,
) -> tuple[TM, TM, TM, dict]:
    """Return ``(S, d/dy S, d/dy' S, trace)`` with a derived remainder."""

    order = c.n
    total = y / y_prime
    derivative_y = y_prime.inverse()
    derivative_y_prime = -(y * y_prime.inverse() * y_prime.inverse())

    current_y = y
    current_y_prime = y_prime
    dy_dy = TM.constant(1, order)
    dy_dyp = TM.zero(order)
    dyp_dy = TM.zero(order)
    dyp_dyp = TM.constant(1, order)

    last_term = total
    last_sensitivity_y = TM.zero(order)
    last_sensitivity_yp = TM.zero(order)
    terms = 1

    c_lower = c.real_lower()
    c_upper = c.bound()
    if not bool(c_lower > 0):
        raise TailError("tail base is not separated from zero")
    if not bool(c_upper <= 2):
        raise TailError("forward invariance needs c <= 2")
    threshold = (2 / c_lower).log()
    # Monotonicity of the tower above the threshold (used by the derivative
    # ratio): c e^u - u is increasing for u >= threshold and equals
    # 2 - threshold there.
    if not bool(threshold < 2):
        raise TailError("tower monotonicity needs log(2/c) < 2")

    closed = None
    for _ in range(1, max_terms):
        inverse_prime = current_y_prime.inverse()
        term = inverse_prime / c
        total = total + term
        sensitivity_y = -(dyp_dy * inverse_prime * inverse_prime / c)
        sensitivity_yp = -(dyp_dyp * inverse_prime * inverse_prime / c)
        derivative_y = derivative_y + sensitivity_y
        derivative_y_prime = derivative_y_prime + sensitivity_yp
        last_term = term
        last_sensitivity_y = sensitivity_y
        last_sensitivity_yp = sensitivity_yp
        terms += 1

        # The ratio to *every* later term is 1/(c y_k) with y_k increasing,
        # so a lower bound on |c y_j| = |c| exp(Re(c y_{j-1})) closes the
        # tail without ever forming the next tower level -- which is what
        # keeps the models representable: one more exponential would widen
        # the base dependence by a factor exp(1e4) over the panel.
        exponent = c * current_y
        exponent_lower = exponent.real_lower()
        tower_lower = c_lower * exponent_lower.exp()
        closable = bool(exponent_lower >= threshold) and bool(tower_lower > 2)
        if closable:
            ratio_upper = 1 / tower_lower
            # The d/dy sensitivities contract by (1 + u_j theta)/(c y_{j+1})
            # <= (1 + u_N)/(c e^{u_N}) at worst, not by 1/(c y_{j+1}); see
            # the module docstring.  Without r_y < 1 the derivative tail is
            # not closed, and neither is the row.
            ratio_y_upper = (1 + exponent_lower) / tower_lower
            closable = bool(ratio_y_upper < 1)
        if closable:
            geometric = ratio_upper / (1 - ratio_upper)
            geometric_y = ratio_y_upper / (1 - ratio_y_upper)
            candidate = last_term.bound() * geometric
            if bool(candidate < stop):
                closed = (
                    ratio_upper,
                    geometric,
                    ratio_y_upper,
                    geometric_y,
                    tower_lower,
                    exponent_lower,
                )
                break

        # One more tower level multiplies the base-panel spread of the model
        # by exp(spread).  Beyond a spread of one the level is no longer
        # representable by any Taylor model of moderate order, so the tail
        # must be closed here -- which is legitimate precisely because the
        # geometric bound above needs only a *scalar* lower bound on the next
        # level, never the level itself.
        spread = exponent.real_upper() - exponent_lower
        if bool(spread > 1):
            if not closable:
                raise TailError(
                    "tower level is not representable and the tail cannot "
                    "be closed at this index"
                )
            closed = (
                ratio_upper,
                geometric,
                ratio_y_upper,
                geometric_y,
                tower_lower,
                exponent_lower,
            )
            break

        # y_{j+1} = exp(c y_j),  y'_{j+1} = c y_{j+1} y'_j, differentiated in
        # the two seeds by the chain rule.
        next_y = (c * current_y).exp()
        next_dy_dy = c * next_y * dy_dy
        next_dy_dyp = c * next_y * dy_dyp
        next_y_prime = c * next_y * current_y_prime
        next_dyp_dy = c * (next_dy_dy * current_y_prime + next_y * dyp_dy)
        next_dyp_dyp = c * (next_dy_dyp * current_y_prime + next_y * dyp_dyp)
        current_y, current_y_prime = next_y, next_y_prime
        dy_dy, dy_dyp = next_dy_dy, next_dy_dyp
        dyp_dy, dyp_dyp = next_dyp_dy, next_dyp_dyp

    if closed is None:
        raise TailError("tail truncation could not be certified")
    (
        ratio_upper,
        geometric,
        ratio_y_upper,
        geometric_y,
        tower_lower,
        exponent_lower,
    ) = closed

    remainder = last_term.bound() * geometric
    remainder_y = last_sensitivity_y.bound() * geometric_y
    remainder_y_prime = last_sensitivity_yp.bound() * geometric

    total = TM(total.c, total.r + remainder, order)
    derivative_y = TM(derivative_y.c, derivative_y.r + remainder_y, order)
    derivative_y_prime = TM(
        derivative_y_prime.c, derivative_y_prime.r + remainder_y_prime, order
    )
    trace = {
        "terms": terms,
        "last_term_upper": last_term.bound(),
        "ratio_upper": ratio_upper,
        "geometric_factor": geometric,
        "ratio_partial_T_upper": ratio_y_upper,
        "geometric_factor_partial_T": geometric_y,
        "tower_lower": tower_lower,
        "invariance_threshold": threshold,
        "remainder": remainder,
        "remainder_partial_T": remainder_y,
        "remainder_partial_T_prime": remainder_y_prime,
    }
    return total, derivative_y, derivative_y_prime, trace
