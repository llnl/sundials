#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
``StabilityFunction``            phi(z) = p(z)/q(z).
``stability_magnitude(phi, Z)``  |phi(Z)|.
``axis_extent(phi, axis)``       extent of the origin-attached stable interval.
``imag_stable_intervals(phi)``   all finite stable intervals on the imaginary axis.
``max_axis_crossing(phi)``       largest |z| where the boundary meets either axis.
``TOL``                          tolerances and resolution floors.
``StabilityWarning``             issued when a result relies on an assumption not confirmed by
                                 the computed data.

The boundary polynomials
------------------------
|phi(z)| <= 1 is equivalent to |p(z)|^2 - |q(z)|^2 <= 0. Restricting z to a half-axis turns the
left-hand side into a real polynomial B in a real parameter:

  * Real axis (z = -t): B(t) = p(-t)^2 - q(-t)^2 = (p(-t) - q(-t)) (p(-t) + q(-t)).

  * Imaginary axis (z = i t): p(i t) and q(i t) split into real parts containing the even powers of
    t and imaginary parts containing the odd powers of t. With y = t^2 this gives the polynomial
    B(y) = Ehat_p(y)^2 + y Ohat_p(y)^2 - Ehat_q(y)^2 - y Ohat_q(y)^2.

The positive roots of B are the axis crossings. Between them, B < 0 means stable and B > 0 means
unstable.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from itertools import pairwise
from typing import Literal, NamedTuple

import numpy as np

from ..utils import EPS, resolved


@dataclass(frozen=True)
class _Tolerances:
    """Tolerances and resolution floors."""

    # Accept a root from np.roots as real despite a small spurious imaginary part.
    root_real_atol: float = 1e-7
    root_real_rtol: float = 1e-6
    # Resolution floor on t; crossings at or below it are treated as the origin.
    root_zero: float = 1e-7
    # Max | |phi| - 1 | for a root of B to count in max_axis_crossing; excludes 0/0 points.
    boundary: float = 1e-6
    # Near-duplicate crossings are merged within this tolerance.
    cluster_rtol: float = 1e-7
    cluster_atol: float = 1e-12


TOL = _Tolerances()


class StabilityWarning(UserWarning):
    """A stability result relies on an assumption not confirmed by the computed data."""


# Extra points sampled inside a segment whose midpoint sign is unresolved.
_SEGMENT_SAMPLES = 7


def _warn(message: str) -> None:
    # stacklevel points at the caller of the public function (axis_extent and friends).
    warnings.warn(message, StabilityWarning, stacklevel=6)


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

    finite_intervals : bounded stable intervals (left, right) in the axis parameter t.
    origin_unbounded : True if the stable region attached to the origin goes to infinity.
    """

    finite_intervals: list[tuple[float, float]]
    origin_unbounded: bool

    @property
    def extent(self) -> float | None:
        """Extent from the origin of the stable interval attached to it.

        None when that interval is unbounded, and 0.0 when the axis leaves the region immediately
        off the origin. The origin-attached interval, when there is one, is the first.
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

    ``coeffs[k]`` multiplies ``x**k`` and ``bounds[k]`` bounds its error. The arithmetic below
    propagates the bounds through the operations that build B. A coefficient is treated as zero
    unless it is resolved against its bound (:func:`~suntools.utils.resolved`).
    """

    coeffs: np.ndarray
    bounds: np.ndarray

    # -- building ------------------------------------------------------------------

    @classmethod
    def zero(cls) -> _BoundedPoly:
        """The zero polynomial. Its sign is unresolved everywhere and it has no roots."""
        return cls(np.zeros(1), np.zeros(1))

    @classmethod
    def from_poly(cls, poly, error) -> _BoundedPoly:
        """From a numpy poly1d and its per-coefficient bounds, both in descending order."""
        return cls(np.asarray(poly.c, float)[::-1], np.asarray(error, float)[::-1])

    # -- arithmetic ----------------------------------------------------------------

    def _combine(self, other: _BoundedPoly, sign: float) -> _BoundedPoly:
        size = max(self.coeffs.size, other.coeffs.size)
        coeffs = _pad(self.coeffs, size) + sign * _pad(other.coeffs, size)
        # Inherited errors add, plus eps for rounding the sum.
        bounds = _pad(self.bounds, size) + _pad(other.bounds, size) + EPS * np.abs(coeffs)
        return _BoundedPoly(coeffs, bounds)

    def __add__(self, other: _BoundedPoly) -> _BoundedPoly:
        return self._combine(other, 1.0)

    def __sub__(self, other: _BoundedPoly) -> _BoundedPoly:
        return self._combine(other, -1.0)

    def __mul__(self, other: _BoundedPoly) -> _BoundedPoly:
        """Product with first-order error |a| eb + ea |b|, plus rounding.

        Each product coefficient sums at most min(len(a), len(b)) terms, so rounding adds at most
        (min + 1) eps times the sum of the term magnitudes.
        """
        magnitude = np.convolve(np.abs(self.coeffs), np.abs(other.coeffs))
        terms = min(self.coeffs.size, other.coeffs.size) + 1
        bounds = np.convolve(np.abs(self.coeffs), other.bounds)
        bounds += np.convolve(self.bounds, np.abs(other.coeffs))
        return _BoundedPoly(
            np.convolve(self.coeffs, other.coeffs), bounds + terms * EPS * magnitude
        )

    def shifted(self) -> _BoundedPoly:
        """x * self; exact, so the bounds shift with the coefficients."""
        return _BoundedPoly(np.pad(self.coeffs, (1, 0)), np.pad(self.bounds, (1, 0)))

    def reflected(self) -> _BoundedPoly:
        """self(-x); exact, so the bounds are unchanged."""
        return _BoundedPoly(self.coeffs * (-1.0) ** np.arange(self.coeffs.size), self.bounds)

    def twisted(self) -> _BoundedPoly:
        """self(i x) with the factor i removed from the odd coefficients.

        i^k = (-1)^(k//2) for even k and i (-1)^(k//2) for odd k, so one real sign pattern covers
        both, and :meth:`alternate` then separates the real (even) and imaginary (odd) parts. Sign
        flips are exact, so the bounds are unchanged.
        """
        twist = (-1.0) ** (np.arange(self.coeffs.size) // 2)
        return _BoundedPoly(self.coeffs * twist, self.bounds)

    def alternate(self, start: int) -> _BoundedPoly:
        """Every other coefficient from *start*, as a polynomial in x^2.

        Returns the zero polynomial if there are none (the odd part of a constant).
        """
        coeffs, bounds = self.coeffs[start::2], self.bounds[start::2]
        return _BoundedPoly(coeffs, bounds) if coeffs.size else _BoundedPoly.zero()

    # -- discarding what is only round-off -----------------------------------------

    def trimmed(self, min_strip: int = 0) -> _BoundedPoly:
        """Drop unresolved coefficients from the top and the order-forced zeros at the bottom.

        Top: a coefficient that cancelled algebraically survives as round-off and would add a root
        near 1/eps. Each coefficient is tested against its *own* bound, not the largest
        coefficient, because genuine leading coefficients can be tiny.

        Bottom: near the origin |phi(i t)|^2 - 1 is O(t^(2m)), so its computed sign is pure
        round-off there. In B this is an exact factor x^m, and dividing it out leaves a constant
        term whose sign is the verdict at the origin (:attr:`_Boundary.origin_stable`).
        ``min_strip`` is the m predicted by the method order. The result spans the resolved
        coefficients at or above ``min_strip``; if there are none, B is zero to within round-off
        and the zero polynomial is returned.
        """
        significant = np.flatnonzero(resolved(self.coeffs, self.bounds))
        significant = significant[significant >= min_strip]
        if significant.size == 0:
            return _BoundedPoly.zero()
        strip, top = int(significant[0]), int(significant[-1]) + 1
        return _BoundedPoly(self.coeffs[strip:top], self.bounds[strip:top])

    # -- interrogation --------------------------------------------------------------

    def value_and_bound(self, x):
        """The polynomial at *x* and a bound on its error; scalar or array.

        The bound is an evaluation at |x| of the coefficient bounds plus the rounding of the
        evaluation itself.
        """
        x = np.asarray(x, dtype=float)
        radius = np.abs(x)
        bound_coeffs = self.bounds + self.coeffs.size * EPS * np.abs(self.coeffs)
        return np.polyval(self.coeffs[::-1], x), np.polyval(bound_coeffs[::-1], radius)

    def signs(self, x):
        """+1 or -1 where the value is resolved against its error bound, 0 elsewhere."""
        value, bound = self.value_and_bound(x)
        usable = np.isfinite(value) & resolved(value, bound)
        return np.where(usable, np.sign(value), 0).astype(int)

    @property
    def leading_sign(self) -> int:
        """Sign of the polynomial for large positive x."""
        return int(np.sign(self.coeffs[-1]))

    def roots(self) -> np.ndarray:
        """Real roots, as candidate crossings.

        np.roots can return inaccurate roots for an ill-conditioned polynomial. A spurious
        root is harmless because verdicts come from signs, not from the root list: an extra
        split between two segments with the same verdict merges away.
        """
        if self.coeffs.size <= 1:
            return np.array([])
        roots = np.roots(self.coeffs[::-1])  # descending, and it drops leading zeros
        real = np.abs(roots.imag) <= TOL.root_real_atol + TOL.root_real_rtol * np.abs(roots.real)
        return roots[real].real


def _imag_mag_sq(bounded: _BoundedPoly) -> _BoundedPoly:
    """|f(i t)|^2 as a polynomial in y = t^2, where f is p or q.

    Writing f(i t) = E(t) + i O(t), E holds the even powers of t and O the odd ones, so
    E(t) = Ehat(t^2), O(t) = t Ohat(t^2), and

        |f(i t)|^2 = Ehat(y)^2 + y Ohat(y)^2.
    """
    split = bounded.twisted()
    even, odd = split.alternate(0), split.alternate(1)
    return even * even + (odd * odd).shifted()


@dataclass(frozen=True)
class _Boundary:
    """The sign of the boundary polynomial B on a half-axis, as a product of factors.

    ``factors`` are polynomials in the boundary coordinate (t or y) with their origin factors
    divided out, so for t > 0 the sign of their product is the sign of B. The axis conventions live
    only in :meth:`point`, :meth:`coordinate`, and :meth:`parameter`.

    B that is zero to within its bounds needs no special case. Its zero factor has unresolved sign
    everywhere, zero leading coefficient, and no roots, which read as "on the boundary", "stable at
    infinity", and "no crossings".
    """

    factors: tuple[_BoundedPoly, ...]
    axis: Literal["real", "imag"]  # validated by _boundary, the only constructor

    def point(self, t):
        """The point z at half-axis parameter t; scalar or array.

        On the real axis z = -t, so t > 0 is the negative real axis. Negative t (the positive
        real axis) is used by :func:`max_axis_crossing` and ignored by the interval classifier.
        """
        t = np.asarray(t, dtype=float)
        return 1j * t if self.axis == "imag" else -t

    def coordinate(self, t):
        """The boundary polynomial's variable at parameter t: y = t^2 (imag) or t (real)."""
        t = np.asarray(t, dtype=float)
        return t * t if self.axis == "imag" else t

    def parameter(self, x) -> np.ndarray:
        """Inverse of :meth:`coordinate`; on the imaginary axis, only y > 0 is kept.

        Callers discard crossings at the origin with TOL.root_zero.
        """
        x = np.asarray(x, dtype=float)
        return np.sqrt(x[x > 0.0]) if self.axis == "imag" else x

    def signs(self, t):
        """+1 outside the region, -1 inside, 0 if any factor is unresolved; array-safe."""
        x = self.coordinate(t)
        product = np.ones(np.shape(x), dtype=int)
        for factor in self.factors:
            product *= factor.signs(x)  # in place, so a 0-d result stays an array
        return product

    def sign(self, t: float) -> int:
        return int(self.signs(t))

    @property
    def origin_stable(self) -> bool:
        """True if the axis is in the region immediately off the origin.

        With the origin factor divided out (B = y^m K on the imaginary axis, B = t K on the
        real axis), this is the sign of K(0). An unresolved sign is on the boundary, which is
        part of the closed region, so it counts as stable.
        """
        return self.sign(0.0) <= 0

    @property
    def far_stable(self) -> bool:
        """True if the axis is in the region far from the origin, from the leading coefficients."""
        return int(np.prod([factor.leading_sign for factor in self.factors])) <= 0

    def roots(self) -> np.ndarray:
        """Candidate crossings from all factors, in the parameter t."""
        return np.concatenate([self.parameter(factor.roots()) for factor in self.factors])


def _check_origin_multiplicity(difference: _BoundedPoly, m: int, order: int | None) -> None:
    """Warn if B(y) has a resolved coefficient below y^m.

    These coefficients are forced to zero based on the method order, so a resolved coefficient
    indicates the order is wrong or the tableau satisfies its order conditions only through severe
    cancellation.
    """
    low = np.flatnonzero(resolved(difference.coeffs[:m], difference.bounds[:m]))
    if low.size:
        k = int(low[0])
        _warn(
            f"order={order} forces B(y) to vanish to order y^{m} on the imaginary axis, but the "
            f"y^{k} coefficient ({difference.coeffs[k]:.3e}) is resolved against its round-off "
            f"bound ({difference.bounds[k]:.1e}); the order metadata or tableau may be wrong"
        )


def _boundary(phi: StabilityFunction, axis: Literal["real", "imag"]) -> _Boundary:
    """Build the boundary for the negative real axis or the imaginary axis.

    The imaginary axis gives one polynomial B(y), with origin multiplicity m set by the method
    order. The real axis keeps B(t) = (p(-t) - q(-t)) (p(-t) + q(-t)) as two factors, which are
    better conditioned to root than their product.
    """
    if axis not in ("real", "imag"):
        raise ValueError(f"axis must be 'real' or 'imag', not {axis!r}")

    p = _BoundedPoly.from_poly(phi.p, phi.p_error)
    q = _BoundedPoly.from_poly(phi.q, phi.q_error)

    if axis == "imag":
        difference = _imag_mag_sq(p) - _imag_mag_sq(q)
        m = phi.origin_multiplicity_y
        _check_origin_multiplicity(difference, m, phi.order)
        return _Boundary((difference.trimmed(m),), axis)

    p_neg, q_neg = p.reflected(), q.reflected()
    factors = ((p_neg - q_neg).trimmed(min_strip=1), (p_neg + q_neg).trimmed())
    return _Boundary(factors, axis)


# ---------------------------------------------------------------------------
# Locating crossings
# ---------------------------------------------------------------------------


def _crossings(boundary: _Boundary) -> list[float]:
    """Crossings on the half-axis in increasing t (the positive roots of B), as floats.

    Crossings come from np.roots and near-duplicate roots are merged, since np.roots splits a
    double root into two copies. Proximity alone cannot separate a split double root from a narrow
    island, so the sign of B at the midpoint decides: unresolved means a double root and is merged,
    resolved means an island and is kept.
    """
    roots = boundary.roots()
    roots = np.sort(roots[roots > TOL.root_zero]).tolist()
    merged: list[float] = []
    for value in roots:
        if merged:
            previous = merged[-1]
            middle = 0.5 * (previous + value)
            close = value - previous <= TOL.cluster_atol + TOL.cluster_rtol * value
            if close and boundary.sign(middle) == 0:
                # The copies split symmetrically, so the midpoint cancels their O(sqrt(eps)) error.
                merged[-1] = middle
                continue
        merged.append(value)
    return merged


# ---------------------------------------------------------------------------
# Classifying stable intervals
# ---------------------------------------------------------------------------


def _segment_stable(boundary: _Boundary, left: float, right: float) -> bool:
    """Stability verdict for the open segment between two consecutive crossings.

    The sign at the midpoint decides. If it is unresolved, more interior points are sampled and the
    first resolved sign decides. A :class:`StabilityWarning` is issued if the resolved samples
    disagree (a missed crossing) or if none resolves, in which case the segment is counted as
    stable.
    """
    mid_sign = boundary.sign(0.5 * (left + right))
    if mid_sign != 0:
        return mid_sign < 0
    samples = np.linspace(left, right, _SEGMENT_SAMPLES + 2)[1:-1]
    signs = boundary.signs(samples)
    resolved_signs = signs[signs != 0]
    if resolved_signs.size == 0:
        _warn(
            f"sign of the boundary polynomial is unresolved throughout the segment "
            f"({left:.6g}, {right:.6g}) on the {boundary.axis} axis; counting it as stable"
        )
        return True
    if np.any(resolved_signs != resolved_signs[0]):
        _warn(
            f"the boundary polynomial changes sign inside the segment ({left:.6g}, "
            f"{right:.6g}) on the {boundary.axis} axis, so a crossing was missed; the "
            f"verdict for this segment is taken from its first resolved sample"
        )
    return bool(resolved_signs[0] < 0)


def _classify_stable_intervals(boundary: _Boundary) -> _AxisIntervals:
    """Group the segments between crossings into maximal stable runs, in the parameter t.

    Every verdict is based on the sign of B. The first segment is judged at the origin
    (:attr:`_Boundary.origin_stable`), the last by :attr:`_Boundary.far_stable`, and the rest by
    :func:`_segment_stable`.
    """
    crossings = _crossings(boundary)
    if not crossings:
        return _AxisIntervals([], boundary.origin_stable)  # uniformly stable or unstable

    # edges[i] opens segment i, and the last segment is unbounded, so there is one verdict per
    # edge.
    edges = [0.0, *crossings]
    verdicts = [boundary.origin_stable]
    verdicts += [_segment_stable(boundary, a, b) for a, b in pairwise(edges[1:])]
    verdicts.append(boundary.far_stable)

    runs: list[tuple[float, float]] = []
    run_start: float | None = None
    # strict=True checks the one-verdict-per-edge invariant.
    for edge, stable in zip(edges, verdicts, strict=True):
        if stable and run_start is None:
            run_start = edge
        elif not stable and run_start is not None:
            runs.append((run_start, edge))
            run_start = None

    # A run still open here is unbounded. If it starts at the origin (edges[0] is exactly 0.0;
    # other edges exceed TOL.root_zero), the whole axis from the origin is stable.
    if run_start == 0.0:
        return _AxisIntervals([], True)
    # Drop detached runs narrower than the resolution floor.
    return _AxisIntervals([(a, b) for a, b in runs if b - a > TOL.root_zero], False)


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def stability_magnitude(phi, Z):
    """|phi(Z)| = |p(Z) / q(Z)| for scalar or array Z."""
    # Silence divide/invalid warnings at poles; inf/nan values are kept for contouring.
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.abs(phi.p(Z) / phi.q(Z))


def _axis_intervals(phi: StabilityFunction, axis: Literal["real", "imag"]) -> _AxisIntervals:
    """Stable intervals on one half-axis; shared by axis_extent and imag_stable_intervals."""
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

    An origin-attached interval has t_left == 0.0; detached islands are included, and
    unbounded stable tails are omitted. Consistent with :func:`axis_extent` by construction.
    """
    return _axis_intervals(phi, "imag").finite_intervals


def max_axis_crossing(phi: StabilityFunction) -> float:
    """Largest |z| at which |phi| = 1 on the real or imaginary axis.

    Includes the positive real axis, since a region can extend into the right half-plane.
    Roots of B where |phi| is not within TOL.boundary of 1 (0/0 points or inaccurate roots)
    are ignored. Used to frame plots.
    """
    largest = 0.0
    for axis in ("real", "imag"):
        boundary = _boundary(phi, axis)
        roots = boundary.roots()
        roots = roots[np.abs(roots) > TOL.root_zero]
        magnitude = stability_magnitude(phi, boundary.point(roots))
        on_boundary = np.isfinite(magnitude) & (np.abs(magnitude - 1.0) <= TOL.boundary)
        if on_boundary.any():
            largest = max(largest, float(np.abs(roots[on_boundary]).max()))
    return largest


__all__ = [
    "StabilityFunction",
    "StabilityWarning",
    "TOL",
    "axis_extent",
    "imag_stable_intervals",
    "max_axis_crossing",
    "stability_magnitude",
]
