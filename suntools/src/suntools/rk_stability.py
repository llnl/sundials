#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# Copyright (c) 2013-2025, Lawrence Livermore National Security
# and Southern Methodist University.
# Copyright (c) 2002-2013, Lawrence Livermore National Security.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------------------

"""Runge-Kutta stability analysis.

Given a method's linear stability function phi(z) = p(z)/q(z), this module analyzes the absolute
stability region { z : |phi(z)| <= 1 }.

Public API
----------
``StabilityFunction``            phi(Z) = p(Z)/q(Z).
``stability_magnitude(phi, Z)``  |phi(Z)|.
``axis_extent(phi, axis)``       extent of the origin-attached stable interval.
``imag_stable_intervals(phi)``   all finite stable intervals on the imaginary axis.
``max_axis_crossing(phi)``       largest |z| where the boundary meets either axis.
``TOL``                          tolerances and resolution floors.

The boundary polynomials
------------------------
|phi(z)| <= 1 is equivalent to |p(z)|^2 - |q(z)|^2 <= 0. Restricting z to a half-axis turns the
left-hand side into a real polynomial B in a real parameter:

  * Real axis (z = -t): B(t) = p(-t)^2 - q(-t)^2 = (p(-t) - q(-t)) (p(-t) + q(-t)).

  * Imaginary axis (z = i t): p(i t) and q(i t) split into real parts containing the even powers of
    t and imaginary parts containing the odd powers of t. With y = t^2 this gives the polynomial
    B(y) = Ehat_p(y)^2 + y Ohat_p(y)^2 - Ehat_q(y)^2 - y Ohat_q(y)^2.

The positive roots of B are the axis crossings, and the sign of B gives the stability of the
segments between the roots. B <= 0 is equivalent to |phi| <= 1, so a segment where B is negative is
stable and one where B is positive is unstable.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import pairwise
from typing import Literal, NamedTuple

import numpy as np

from .utils import EPS, resolved


@dataclass(frozen=True)
class _Tolerances:
    """Tolerances and resolution floors."""

    # np.roots works in the complex plane even for real polynomials, so a real root could be
    # returned with a small spurious imaginary part.
    root_real_atol: float = 1e-7
    root_real_rtol: float = 1e-6
    # Resolution floor on t, and the discard threshold for the trivial crossing at t=0. It is
    # applied in t by every consumer of a root list, so it needs no y = t^2 counterpart.
    root_zero: float = 1e-7
    # How close |phi| must be to 1 for a root of B to be reported as a crossing of the plot's
    # framing measure; its only job is to exclude poles and 0/0 points.
    boundary: float = 1e-6
    # Near-duplicate crossings are merged within this tolerance.
    cluster_rtol: float = 1e-7
    cluster_atol: float = 1e-12


TOL = _Tolerances()


@dataclass(frozen=True)
class StabilityFunction:
    """phi(z) = p(z) / q(z)"""

    p: np.poly1d
    q: np.poly1d
    p_error: np.ndarray  # round-off bounds on p
    q_error: np.ndarray  # round-off bounds on q
    order: int | None = None  # RK method order

    @property
    def origin_multiplicity_y(self) -> int:
        """Multiplicity m of y = 0 as a root of the boundary polynomial B(y), y = t^2.

        A method of order r has phi(z) = exp(z) + C z^(r+1) + ..., so |phi(i t)|^2 - 1 has a zero
        of order at least r+1 at t = 0, as does B = |q(i t)|^2 (|phi(i t)|^2 - 1) since q(0) = 1. B
        is even in t, so that order is 2m >= r+1, giving the floor m >= ceil((r+1)/2).
        """
        return 0 if self.order is None else (self.order + 2) // 2


class _AxisIntervals(NamedTuple):
    """Stability classification along one half-axis.

    finite_intervals : bounded stable intervals (left, right) in the axis parameter t. A
    plain list, and not defended against mutation: :func:`_axis_intervals` builds a fresh
    classification on every call, so no two readers ever hold the same one. Caching them
    would have to revisit that.

    origin_unbounded : True if the stable region attached to the origin goes to infinity.

    Both public readers take what they report from one of these, so the interval list and
    the reported extent are two views of a single classification and cannot disagree.
    """

    finite_intervals: list[tuple[float, float]]
    origin_unbounded: bool

    @property
    def extent(self) -> float | None:
        """Extent from the origin of the stable interval attached to it.

        None when that interval is unbounded, and 0.0 when the axis leaves the region
        immediately off the origin. The origin-attached interval, when there is one, is the
        first and starts at exactly 0.0 -- :func:`_classify_stable_intervals` opens it on the
        origin edge -- so the test is an identity rather than a tolerance.
        """
        if self.origin_unbounded:
            return None
        if self.finite_intervals:
            left, right = self.finite_intervals[0]
            if left == 0.0:
                return right
        return 0.0


# ---------------------------------------------------------------------------
# Polynomials that carry their own round-off
# ---------------------------------------------------------------------------


def _pad(values, size):
    return np.pad(np.asarray(values, float), (0, size - np.size(values)))


@dataclass(frozen=True)
class _BoundedPoly:
    """A real polynomial that carries a bound on the round-off in each coefficient.

    ``coeffs[k]`` multiplies ``x**k`` and ``bounds[k]`` bounds the error already present in it.
    The arithmetic below propagates them through the squarings and differences that build B. A
    coefficient counts as zero when it falls below its own bound; that test is
    :func:`~suntools.utils.resolved`, shared with the code that builds p and q so the two cannot
    disagree.

    Coefficients ascend because the origin end -- where the order conditions force the
    cancellation that matters -- is then simply the front.
    """

    coeffs: np.ndarray
    bounds: np.ndarray

    # -- building ------------------------------------------------------------------

    @classmethod
    def zero(cls, bound: float = 0.0) -> _BoundedPoly:
        """The zero polynomial, carrying *bound* as the round-off in its one coefficient.

        The default is the exact zero that :meth:`alternate` wants: a missing parity is
        missing exactly, not to within some error. :meth:`trimmed` passes its own bound
        instead, since a polynomial that is zero only to within its round-off still knows
        how large that round-off was.
        """
        return cls(np.zeros(1), np.full(1, bound))

    @classmethod
    def from_poly(cls, poly, error) -> _BoundedPoly:
        """From a numpy poly1d and the matching descending per-coefficient bounds."""
        return cls(np.asarray(poly.c, float)[::-1], np.asarray(error, float)[::-1])

    # -- arithmetic ----------------------------------------------------------------

    def _combine(self, other: _BoundedPoly, sign: float) -> _BoundedPoly:
        size = max(self.coeffs.size, other.coeffs.size)
        coeffs = _pad(self.coeffs, size) + sign * _pad(other.coeffs, size)
        # Inherited error adds; the addition itself rounds to within eps of its own result.
        bounds = _pad(self.bounds, size) + _pad(other.bounds, size) + EPS * np.abs(coeffs)
        return _BoundedPoly(coeffs, bounds)

    def __add__(self, other: _BoundedPoly) -> _BoundedPoly:
        return self._combine(other, 1.0)

    def __sub__(self, other: _BoundedPoly) -> _BoundedPoly:
        return self._combine(other, -1.0)

    def __mul__(self, other: _BoundedPoly) -> _BoundedPoly:
        """A coefficient of the product is a sum of at most min(m, n) terms a_i b_j, so the
        first-order error is |a| eb + ea |b| and the arithmetic adds at most (min + 1) eps
        times the sum of the term magnitudes.
        """
        magnitude = np.convolve(np.abs(self.coeffs), np.abs(other.coeffs))
        terms = min(self.coeffs.size, other.coeffs.size) + 1
        bounds = np.convolve(np.abs(self.coeffs), other.bounds)
        bounds += np.convolve(self.bounds, np.abs(other.coeffs))
        return _BoundedPoly(
            np.convolve(self.coeffs, other.coeffs), bounds + terms * EPS * magnitude
        )

    def shifted(self) -> _BoundedPoly:
        """x * self: exact in floating point, so the bounds only move up with it."""
        return _BoundedPoly(np.pad(self.coeffs, (1, 0)), np.pad(self.bounds, (1, 0)))

    def reflected(self) -> _BoundedPoly:
        """self(-x): negating alternate coefficients is exact, so the bounds pass through."""
        return _BoundedPoly(self.coeffs * (-1.0) ** np.arange(self.coeffs.size), self.bounds)

    def twisted(self) -> _BoundedPoly:
        """self(i x) with the i factored out of the odd coefficients.

        i^k is (-1)^(k//2) for even k and i (-1)^(k//2) for odd k, so one real sign pattern
        serves both parities and :meth:`alternate` can then take the two of them apart as the
        real and imaginary parts. Like :meth:`reflected` this is a sign flip and nothing more,
        so the bounds pass through untouched.
        """
        twist = (-1.0) ** (np.arange(self.coeffs.size) // 2)
        return _BoundedPoly(self.coeffs * twist, self.bounds)

    def alternate(self, start: int) -> _BoundedPoly:
        """Every other coefficient from *start*, as a polynomial in x^2.

        A constant polynomial has no odd part, and the zero polynomial is what that means.
        """
        coeffs, bounds = self.coeffs[start::2], self.bounds[start::2]
        return _BoundedPoly(coeffs, bounds) if coeffs.size else _BoundedPoly.zero()

    # -- discarding what is only round-off -----------------------------------------

    def trimmed(self, min_strip: int = 0) -> _BoundedPoly:
        """Drop both ends that are round-off rather than signal.

        The two ends are not symmetric and both matter.

        At the top, a coefficient that cancelled algebraically survives as round-off and
        would contribute a root of order 1/eps. Testing it against its *own* bound rather
        than against the polynomial's largest coefficient is what makes that safe: the
        genuine leading coefficient of the boundary polynomial of an s-stage explicit
        method is (1/s!)^2, which for s >= 11 is already below eps times an O(1) constant
        term, and any scale-relative test would throw it away.

        At the origin, the x^m that the order conditions guarantee has to come out, or the
        sign near the origin is undecidable. Near the origin on the imaginary axis
        |phi(i t)| - 1 is O(t^(2m)) with 2m >= order + 1, so for any method of order >= 2
        there is a neighbourhood in which the computed sign of |phi| - 1 is pure round-off
        -- no tolerance can rescue it, because the information is not in the double. In B
        that cancellation is an exact factor rather than a near-cancellation, so dividing
        it out recovers the verdict exactly (:attr:`_Boundary.origin_stable` reads it).
        ``min_strip`` is the multiplicity theory predicts, and the larger of it and the
        measured one wins: for a well-conditioned tableau the bounds find m on their own,
        while for one whose order conditions hold only through cancellation the round-off
        at the origin can swamp the test and leave genuine-looking coefficients that theory
        says are exactly zero.

        Everything downstream assumes it is looking at a polynomial in this form. When
        nothing anywhere is resolved the result is the zero polynomial; what that means is
        the caller's business (see :class:`_Boundary`).
        """
        significant = np.flatnonzero(resolved(self.coeffs, self.bounds))
        if significant.size == 0:
            return _BoundedPoly.zero(float(self.bounds[0]))
        top = int(significant[-1]) + 1
        strip = min(max(int(significant[0]), min_strip), top - 1)
        return _BoundedPoly(self.coeffs[strip:top], self.bounds[strip:top])

    # -- interrogation --------------------------------------------------------------

    def value_and_bound(self, x):
        """The polynomial at *x* and a bound on the error in that value; scalar or array.

        Both are Horner evaluations, the bound at |x| so that it is an upper bound term by
        term. The bound's coefficients are the round-off each one arrived with plus the
        round-off of summing this particular series; the two are folded together first
        because the bound is linear in them, which leaves one Horner pass rather than two.
        """
        x = np.asarray(x, dtype=float)
        radius = np.abs(x)
        bound_coeffs = self.bounds + self.coeffs.size * EPS * np.abs(self.coeffs)
        return np.polyval(self.coeffs[::-1], x), np.polyval(bound_coeffs[::-1], radius)

    def signs(self, x):
        """+1, -1, or 0 where the value does not stand clear of its own round-off.

        The evaluation carries its bound along, so a sign is reported only where it is
        resolved. Where it is not, x is within round-off of a root, and the verdict there
        cannot move an interval endpoint anyway.
        """
        value, bound = self.value_and_bound(x)
        usable = np.isfinite(value) & resolved(value, bound)
        return np.where(usable, np.sign(value), 0).astype(int)

    @property
    def leading_sign(self) -> int:
        """Sign of the polynomial for large positive x."""
        return int(np.sign(self.coeffs[-1]))

    def roots(self) -> np.ndarray:
        """Real roots.

        Roots are candidates only: np.roots on an ill-conditioned high-degree polynomial
        can return values that do not satisfy the equation. A stray one is harmless,
        because the stability verdict does not come from the root list -- an extra
        partition point between two points of the same sign simply merges away.
        """
        if self.coeffs.size <= 1:
            return np.array([])
        roots = np.roots(self.coeffs[::-1])  # descending, and it drops leading zeros
        real = np.abs(roots.imag) <= TOL.root_real_atol + TOL.root_real_rtol * np.abs(roots.real)
        return roots[real].real


def _imag_mag_sq(bounded: _BoundedPoly) -> _BoundedPoly:
    """|f(i t)|^2 as a polynomial in y = t^2, where f is p or q.

    Split f(i t) into its real and imaginary parts, f(i t) = E(t) + i O(t) with E
    and O real; then |f(i t)|^2 = E(t)^2 + O(t)^2. E collects the even powers of t and O
    the odd ones (:meth:`_BoundedPoly.twisted` applies the i^k sign pattern that makes this
    a split by parity), so E(t) = Ehat(t^2) and O(t) = t Ohat(t^2) and

        |f(i t)|^2 = Ehat(y)^2 + y Ohat(y)^2,    y = t^2.

    These are the module docstring's Ehat_p, Ohat_p for f = p and Ehat_q, Ohat_q for f = q;
    B(y) is the difference of the two results. Working in y halves the degree handed to the
    root finder.
    """
    split = bounded.twisted()
    even, odd = split.alternate(0), split.alternate(1)
    return even * even + (odd * odd).shifted()


@dataclass(frozen=True)
class _Boundary:
    """The sign of the boundary polynomial B on a half-axis, as a product of factors.

    ``factors`` are polynomials in the boundary coordinate (t or y), each with its origin
    factor already divided out, so for t > 0 the sign of the product is the sign of B.
    ``axis`` fixes the half-axis, and with it the three mappings below --
    :meth:`point`, :meth:`coordinate`, and its inverse :meth:`parameter`. Those are the
    only statements of the convention anywhere in the module; nothing else branches on the
    axis, so the two halves of it cannot drift apart.

    A boundary polynomial that is zero to within its own bounds needs no special case: it
    arrives here as a zero factor, whose sign is unresolved everywhere, whose leading
    coefficient is zero, and which has no roots. The three fall out as "on the boundary
    everywhere", "stable at infinity", and "no crossings" -- which is exactly right, and
    is what a symmetric method looks like. Implicit midpoint gives a difference that is
    identically zero, the trapezoidal rule one of 1.1e-16 against terms of order 1.
    """

    factors: tuple[_BoundedPoly, ...]
    # :func:`_boundary` is the only constructor, and it validates the axis it was handed
    # before doing any work, so there is nothing left for a __post_init__ to catch.
    axis: Literal["real", "imag"]

    def point(self, t):
        """The point z in the complex plane at half-axis parameter t; scalar or array.

        The real axis runs z = -t, so that t increases into the left half-plane where the
        stable interval is. A *negative* t is therefore a point on the positive real axis,
        which the plot's framing measure wants and the interval classifier discards.
        """
        t = np.asarray(t, dtype=float)
        return 1j * t if self.axis == "imag" else -t

    def coordinate(self, t):
        """The boundary polynomial's variable at half-axis parameter t.

        The imaginary axis works in y = t^2, which halves the degree handed to the root
        finder; the real axis works in t directly.
        """
        t = np.asarray(t, dtype=float)
        return t * t if self.axis == "imag" else t

    def parameter(self, x) -> np.ndarray:
        """The half-axis parameters whose coordinate is *x*: the inverse of :meth:`coordinate`.

        On the imaginary axis y = t^2 has a real t only for positive y, which is all the
        filter here enforces. Dropping the crossing that sits at the origin is left to the
        consumers, which apply TOL.root_zero in t: a y that survives this filter but falls
        below root_zero^2 maps to a t they discard anyway.
        """
        x = np.asarray(x, dtype=float)
        return np.sqrt(x[x > 0.0]) if self.axis == "imag" else x

    def signs(self, t):
        """+1 outside the region, -1 inside, 0 where the sign is not resolved; array-safe.

        A single unresolved factor makes the product unresolved, whatever the others say.
        """
        x = self.coordinate(t)
        product = np.ones(np.shape(x), dtype=int)
        for factor in self.factors:
            product *= factor.signs(x)  # in place, so a 0-d result stays an array
        return product

    def sign(self, t: float) -> int:
        return int(self.signs(t))

    def stable(self, t: float) -> bool:
        """Is the half-axis point at parameter t in the closed stability region?

        An unresolved sign means the point is within round-off of the boundary, which is
        part of the closed region, so it counts as stable.
        """
        return self.sign(t) <= 0

    @property
    def origin_stable(self) -> bool:
        """Is the axis in the region immediately off the origin?

        Exact, not sampled. Each factor arrives with its origin power already divided out,
        leaving K -- B = y^m K on the imaginary axis, and B = t K on the real axis, where
        phi(0) = 1 makes t = 0 a simple root -- so this is the sign of a product of
        constant terms K(0), each an ordinary O(1) number. For Heun's method B(y) = y^2 / 4,
        so K = 1/4 > 0 and the method is outside the region for every t > 0; for classical
        RK4 B(y) = y^3 (y - 8) / 576, so K(0) < 0 and the crossing is at y = 8, t = 2 sqrt 2.
        """
        return self.stable(0.0)

    @property
    def far_sign(self) -> int:
        """The sign arbitrarily far out along the axis.

        The sign of a polynomial at infinity is the sign of its leading coefficient, so
        there is no need to evaluate a high-degree polynomial at a large argument.
        """
        return int(np.prod([factor.leading_sign for factor in self.factors]))

    @property
    def far_stable(self) -> bool:
        """Is the axis in the region arbitrarily far out?"""
        return self.far_sign <= 0

    def roots(self) -> np.ndarray:
        """Every real candidate crossing, in the axis parameter t.

        The factors live in the boundary coordinate; :meth:`parameter` brings them back to
        t, so everything downstream works in the one parameter. Non-finite values are
        filtered here, which is what lets the consumers take finiteness for granted.
        """
        # Both branches of :func:`_boundary` supply at least one factor, so an empty
        # ``factors`` is not a state this class has to represent; a zero factor is (see the
        # class docstring), and that is a different thing.
        roots = np.concatenate([self.parameter(factor.roots()) for factor in self.factors])
        return roots[np.isfinite(roots)]


def _boundary(phi: StabilityFunction, axis: Literal["real", "imag"]) -> _Boundary:
    """Build the boundary object for the negative real axis or the imaginary axis.

    The imaginary axis gives one polynomial, B(y) in y = t^2, whose origin multiplicity m
    the method's order predicts. The real axis gives two polynomials in t, because
    B(t) = p(-t)^2 - q(-t)^2 factors exactly as (p(-t) - q(-t)) (p(-t) + q(-t)) and rooting
    the two degree-<=s factors is better conditioned than rooting their degree-2s product.

    The axis is checked here, before any of that, because this is the only place a
    :class:`_Boundary` is built and the only place an axis coming from the public API is
    still just an argument rather than a branch someone has already fallen through.
    """
    if axis not in ("real", "imag"):
        raise ValueError(f"axis must be 'real' or 'imag', not {axis!r}")

    p = _BoundedPoly.from_poly(phi.p, phi.p_error)
    q = _BoundedPoly.from_poly(phi.q, phi.q_error)

    if axis == "imag":
        difference = _imag_mag_sq(p) - _imag_mag_sq(q)
        return _Boundary((difference.trimmed(phi.origin_multiplicity_y),), axis)

    p_neg, q_neg = p.reflected(), q.reflected()
    # phi(0) = 1 makes t = 0 a root of p(-t) - q(-t); order >= 1 makes it a simple one,
    # since phi(-t) - 1 = -t + O(t^2). The other factor is 2 at the origin, so there is
    # nothing to strip from it.
    factors = ((p_neg - q_neg).trimmed(min_strip=1), (p_neg + q_neg).trimmed())
    return _Boundary(factors, axis)


# ---------------------------------------------------------------------------
# Locating crossings
# ---------------------------------------------------------------------------


def _crossings(boundary: _Boundary) -> list[float]:
    """Every crossing on the half-axis, in the parameter t: the positive roots of B.

    Returned as ordinary floats, in increasing order. The merge below has to build them one
    at a time anyway, since each comparison is against the running merged value rather than
    against the raw root, so there is nothing for an array to vectorize and the caller wants
    floats in the end.

    There is deliberately no second, independent source of crossings and no bisection
    refinement of these. A sweep for sign changes can only ever find odd-multiplicity
    roots, a subset of what the root finder returns, so it cannot cover a root-finder
    failure that matters; and the root finder does not fail here. Across the ARKODE
    tables, Taylor truncations to s = 30 (where B's coefficients span 1e65), and several
    hundred random explicit and implicit methods, a sweep contributed no crossing np.roots
    had not already found, and bisecting B around its roots moved no endpoint by more than
    round-off.

    Near-duplicates are merged: a double root arrives from np.roots as two split copies,
    and the two factors of the real-axis polynomial can report the same crossing. But
    proximity alone cannot tell a split double root from a genuinely narrow island -- a
    random degree-8 polynomial readily produces one a few 1e-6 wide out at t = 70, well
    inside any sensible relative tolerance there. The boundary polynomial settles it, as
    it settles everything else here: between a split double root B is within round-off of
    zero and its sign is unresolved, while across a narrow island the sign is resolved and
    opposite. Only the former is merged.
    """
    roots = boundary.roots()  # already finite: _Boundary.roots filters
    roots = np.sort(roots[roots > TOL.root_zero])
    if roots.size == 0:
        return []
    merged = [float(roots[0])]
    for value in roots[1:]:
        previous = merged[-1]
        threshold = TOL.cluster_atol + TOL.cluster_rtol * max(abs(value), abs(previous))
        middle = 0.5 * (previous + float(value))
        if abs(value - previous) <= threshold and boundary.sign(middle) == 0:
            # A double root splits symmetrically about the true value, so the midpoint
            # cancels the O(sqrt(eps)) error that either copy alone carries.
            merged[-1] = middle
        else:
            merged.append(float(value))
    return merged


# ---------------------------------------------------------------------------
# Classifying stable intervals
# ---------------------------------------------------------------------------


def _classify_stable_intervals(boundary: _Boundary) -> _AxisIntervals:
    """Turn crossings into maximal stable runs along a half-axis, in the parameter t.

    The crossings cut the half-axis into segments, each with a single verdict, and the
    stable ones are then grouped into runs. Every verdict is the sign of the boundary
    polynomial; none is a comparison of |phi| against 1.

    The segment touching the origin is judged at the origin itself rather than at its
    midpoint. That is where K(0) sits (:attr:`_Boundary.origin_stable`), and it is the one
    place the verdict is exact -- a midpoint may still be inside a high-order tangency the
    double cannot resolve.
    """
    crossings = _crossings(boundary)
    if not crossings:
        return _AxisIntervals([], boundary.origin_stable)  # uniformly stable or unstable

    # edges[i] opens segment i; the last segment is the unbounded tail past the last
    # crossing, so there are exactly as many verdicts as edges. The pairing below walks
    # successive interior crossings, so its two ends differ in length by one on purpose --
    # pairwise says that, where a zip of two offset slices only implied it.
    edges = [0.0, *crossings]
    verdicts = [boundary.origin_stable]
    verdicts += [boundary.stable(0.5 * (a + b)) for a, b in pairwise(edges[1:])]
    verdicts.append(boundary.far_stable)

    runs: list[tuple[float, float]] = []
    run_start: float | None = None
    # strict=True turns the count claim above into a check: the three lines that build
    # verdicts have to keep matching edges, and a refactor that broke that would otherwise
    # silently drop the verdict on the final segment.
    for edge, stable in zip(edges, verdicts, strict=True):
        if stable and run_start is None:
            run_start = edge
        elif not stable and run_start is not None:
            runs.append((run_start, edge))
            run_start = None

    # A run still open when the loop ends is one that reaches infinity, so it has no finite
    # right endpoint to report and never enters ``runs``. All that is left to ask is whether
    # it is the run attached to the origin. That is an identity, not a tolerance: edges[0] is
    # the origin exactly, and every other edge is a crossing already filtered to
    # > TOL.root_zero by :func:`_crossings`.
    if run_start == 0.0:  # None, meaning no unbounded run, compares false as it should
        return _AxisIntervals([], True)
    # The width filter is the resolution floor applied to a run rather than to a crossing. It
    # can only ever touch a detached run: an origin-attached one runs from 0.0 to a crossing
    # that _crossings has already filtered to > TOL.root_zero. What it drops is a pair of
    # crossings too close together to report as an interval but far enough apart that the
    # merge in _crossings kept them separate, having resolved an opposite sign between them.
    # A genuine narrow island is wider than this: the 1e-6-wide one that _crossings discusses
    # clears the floor by an order of magnitude.
    return _AxisIntervals([(a, b) for a, b in runs if b - a > TOL.root_zero], False)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def stability_magnitude(phi, Z):
    """|phi(Z)| = |p(Z) / q(Z)| for scalar or array Z."""
    # Divide and invalid warnings are suppressed to keep batch plotting quiet when poles are
    # present; the values (including inf/nan) are preserved for plotting contours.
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.abs(phi.p(Z) / phi.q(Z))


def _axis_intervals(phi: StabilityFunction, axis: Literal["real", "imag"]) -> _AxisIntervals:
    """The half-axis classification that the two readers below are views of."""
    return _classify_stable_intervals(_boundary(phi, axis))


def axis_extent(phi: StabilityFunction, axis: Literal["real", "imag"]) -> float | None:
    """Extent from the origin of the stable interval attached to it.

    ``axis`` is 'real' (z = -t) or 'imag' (z = i t). Returns the right endpoint t, or
    None if that interval is unbounded (e.g. A-stable along the negative real axis), or
    0.0 if the axis leaves the region immediately off the origin.
    """
    return _axis_intervals(phi, axis).extent


def imag_stable_intervals(phi: StabilityFunction) -> list[tuple[float, float]]:
    """Finite stable intervals (t_left, t_right) on the positive imaginary axis.

    Intervals attached to the origin appear with t_left == 0.0; detached islands are
    included too. Unbounded stable tails are omitted: they have no endpoint to annotate.

    This and :func:`axis_extent` are two views of one :class:`_AxisIntervals`, so the
    interval list and the reported extent cannot disagree.
    """
    return list(_axis_intervals(phi, "imag").finite_intervals)


def max_axis_crossing(phi: StabilityFunction) -> float:
    """Largest |z| at which |phi| = 1 on the real or imaginary axis.

    Crossings on the *positive* real axis count too. A stability region can reach into the
    right half-plane, and a plot framed only on the left would clip it.

    A root of the boundary polynomial has |p| = |q|, which is |phi| = 1 unless both vanish.
    The magnitude test below keeps only the roots where it really is 1, which drops those
    removable singularities along with any root the root finder placed too far from its true
    location to frame a plot around. It runs on the whole root array at once: every step
    from :meth:`_Boundary.point` through :func:`stability_magnitude` is array-safe.
    """
    largest = 0.0
    for axis in ("real", "imag"):
        boundary = _boundary(phi, axis)
        roots = boundary.roots()
        roots = roots[np.abs(roots) > TOL.root_zero]
        if roots.size == 0:
            continue
        magnitude = stability_magnitude(phi, boundary.point(roots))
        on_boundary = np.isfinite(magnitude) & (np.abs(magnitude - 1.0) <= TOL.boundary)
        if on_boundary.any():
            largest = max(largest, float(np.abs(roots[on_boundary]).max()))
    return largest


__all__ = [
    "StabilityFunction",
    "TOL",
    "axis_extent",
    "imag_stable_intervals",
    "max_axis_crossing",
    "stability_magnitude",
]
