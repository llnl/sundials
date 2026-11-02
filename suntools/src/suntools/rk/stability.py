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
``r_infinity(phi)``              |phi(infinity)|.
``is_a_stable(phi)``             A-stability, verified from the imaginary axis and the poles.
``is_l_stable(phi)``             A-stable and |phi(infinity)| = 0.
``a_alpha_angle(phi)``           the A(alpha) angle in degrees.
``threshold_factor(phi)``        linear threshold factor (polynomial phi only).
``order_star_function(phi, Z)``  log(|phi(Z)| / |exp(Z)|).
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

Stiff properties
----------------
A-stability needs the imaginary axis in the region and no pole in the open left half-plane. A(alpha)
applies the same machinery to the ray z = -r e^{-i theta}, giving a per-angle verdict, and bisects
theta: a bounded unstable component contains a pole, so the stable angles form one interval [0,
alpha] once the guard scan has excluded detached sectors at infinity. |phi(infinity)| follows from
the degrees of the resolved coefficients.
"""

from __future__ import annotations

import math
import os
import warnings
from dataclasses import dataclass, field
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
    # Looser merge for an np.roots split of a double root (about sqrt(eps * condition) relative),
    # used only when the whole gap is unresolved.
    cluster_unresolved_rtol: float = 1e-5
    # Resolution of the A(alpha) angle search, in degrees.
    angle: float = 1e-9
    # A pole with |Re| <= pole_rtol |pole| counts as on the imaginary axis.
    pole_rtol: float = 1e-6
    # Relative resolution of the radius searches (SSP coefficient, threshold factor).
    radius_rtol: float = 1e-12


TOL = _Tolerances()


class StabilityWarning(UserWarning):
    """A stability result relies on an assumption not confirmed by the computed data."""


# Extra points sampled inside a segment whose midpoint sign is unresolved.
_SEGMENT_SAMPLES = 7
# Multiples of the last crossing sampled beyond it: 1.125x to 9x.
_FAR_SAMPLE_FACTORS = 1.0 + 2.0 ** np.arange(-3, 4)
# Angular step, in degrees, of the guard scan that precedes the A(alpha) bisection.
_A_ALPHA_SCAN_STEP_DEG = 1.0
# Largest radius searched for the SSP coefficient before reporting it as unbounded.
_RADIUS_SEARCH_MAX = 2**20


# Warnings skip frames from the suntools package (this file's parent's parent directory).
_WARN_SKIP_PREFIXES = (os.path.dirname(os.path.dirname(__file__)) + os.sep,)


def _warn(message: str) -> None:
    # Warnings point at the first frame outside suntools.
    warnings.warn(message, StabilityWarning, skip_file_prefixes=_WARN_SKIP_PREFIXES)


@dataclass(frozen=True)
class StabilityFunction:
    """phi(z) = p(z) / q(z)

    p_error and q_error are optional per-coefficient round-off bounds on p and q. They are in
    descending order, have the same length as p.c and q.c, and default to EPS*|c|. The boundary
    polynomials are cached on first use, so the coefficients must not be mutated in place.
    """

    p: np.poly1d
    q: np.poly1d
    order: int | None = None  # RK method order
    # None (the default) means EPS*|c|; __post_init__ always stores an ndarray.
    p_error: np.ndarray | None = field(default=None, kw_only=True, compare=False)
    q_error: np.ndarray | None = field(default=None, kw_only=True, compare=False)
    _boundaries: dict[str, _Boundary] = field(
        default_factory=dict, init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        for name, poly in (("p", self.p), ("q", self.q)):
            value = getattr(self, f"{name}_error")
            if value is None:
                error = EPS * np.abs(poly.c)
            else:
                error = np.asarray(value, float)
            if error.shape != poly.c.shape:
                raise ValueError(
                    f"{name}_error must have the same shape as {name}.c: got length "
                    f"{error.size} for {len(poly.c)} coefficients"
                )
            object.__setattr__(self, f"{name}_error", error)

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


@dataclass(frozen=True, eq=False)
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


@dataclass(frozen=True, eq=False)
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
    # "real" and "imag" are validated by _boundary before it calls _build_boundary; "ray" is
    # internal, for A(alpha).
    axis: Literal["real", "imag", "ray"]
    # Ray only: the angle, in radians, between the ray and the negative real axis (upper half).
    angle: float = 0.0

    @property
    def label(self) -> str:
        """The boundary's name in warning messages."""
        if self.axis == "ray":
            return f"ray at {math.degrees(self.angle):.6g} deg"
        return f"{self.axis} axis"

    def point(self, t):
        """The point z at half-axis parameter t; scalar or array.

        On the real axis z = -t, so t > 0 is the negative real axis. Negative t (the positive
        real axis) is used by :func:`max_axis_crossing` and ignored by the interval classifier.
        """
        t = np.asarray(t, dtype=float)
        if self.axis == "ray":
            return -t * np.exp(-1j * self.angle)
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
    """Return the cached boundary for the negative real axis or the imaginary axis."""
    if axis not in ("real", "imag"):
        raise ValueError(f"axis must be 'real' or 'imag', not {axis!r}")

    if axis not in phi._boundaries:
        phi._boundaries[axis] = _build_boundary(phi, axis)
    return phi._boundaries[axis]


def _build_boundary(phi: StabilityFunction, axis: Literal["real", "imag"]) -> _Boundary:
    """Build the boundary for the negative real axis or the imaginary axis.

    The imaginary axis gives one polynomial B(y), with origin multiplicity m set by the method
    order. The real axis keeps B(t) = (p(-t) - q(-t)) (p(-t) + q(-t)) as two factors, which are
    better conditioned to root than their product.
    """
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


def _ray_mag_sq(f: _BoundedPoly, theta: float) -> _BoundedPoly:
    """|f(z)|^2 as a polynomial in r on the ray z = -r e^{-i theta}, where f is p or q.

    With f(z) = sum_k f_k z^k, the real part has coefficients f_k (-1)^k cos(k theta) and the
    imaginary part -f_k (-1)^k sin(k theta). The angle enters the bounds as an absolute error of
    (k + 1) eps in each cosine and sine.
    """
    k = np.arange(f.coeffs.size)
    base = f.coeffs * (-1.0) ** k
    slack = (k + 2) * EPS * np.abs(f.coeffs)
    cos, sin = np.cos(k * theta), np.sin(k * theta)
    real = _BoundedPoly(base * cos, f.bounds * np.abs(cos) + slack)
    imag = _BoundedPoly(-base * sin, f.bounds * np.abs(sin) + slack)
    return real * real + imag * imag


def _build_ray_boundary(phi: StabilityFunction, angle_deg: float) -> _Boundary:
    """Boundary on the ray at *angle_deg* above the negative real axis, as one factor.

    G(r) = |p(z)|^2 - |q(z)|^2 has the origin factor r divided out, as on the real axis. The result
    is not cached because the angle varies.
    """
    theta = math.radians(angle_deg)
    p = _BoundedPoly.from_poly(phi.p, phi.p_error)
    q = _BoundedPoly.from_poly(phi.q, phi.q_error)
    difference = _ray_mag_sq(p, theta) - _ray_mag_sq(q, theta)
    return _Boundary((difference.trimmed(min_strip=1),), "ray", theta)


def _ray_stable(phi: StabilityFunction, angle_deg: float) -> bool:
    """True if the whole ray from the origin at *angle_deg* lies in the region."""
    return _classify_stable_intervals(_build_ray_boundary(phi, angle_deg)).origin_unbounded


# ---------------------------------------------------------------------------
# Locating crossings
# ---------------------------------------------------------------------------


def _interior_signs(boundary: _Boundary, left: float, right: float) -> np.ndarray:
    """Signs of B at the interior sample points of (left, right); the 4th is the midpoint."""
    return boundary.signs(np.linspace(left, right, _SEGMENT_SAMPLES + 2)[1:-1])


def _crossings(boundary: _Boundary) -> list[float]:
    """Crossings on the half-axis in increasing t (the positive roots of B), as floats.

    Crossings come from np.roots and near-duplicate roots are merged, since np.roots splits a
    double root into two copies. Proximity alone cannot separate a split double root from a narrow
    island, so the sign of B at the midpoint decides: unresolved means a double root and is merged,
    resolved means an island and is kept. A wider gap (up to ``TOL.cluster_unresolved_rtol``) is
    also merged when B is unresolved throughout it, because a double root split by np.roots looks
    like that. A gap with any resolved sample is never merged.
    """
    roots = boundary.roots()
    roots = np.sort(roots[roots > TOL.root_zero]).tolist()
    merged: list[float] = []
    for value in roots:
        if merged:
            previous = merged[-1]
            middle = 0.5 * (previous + value)
            gap = value - previous
            tight = gap <= TOL.cluster_atol + TOL.cluster_rtol * value
            loose = gap <= TOL.cluster_unresolved_rtol * value
            if (
                (tight or loose)
                and boundary.sign(middle) == 0
                and (tight or not _interior_signs(boundary, previous, value).any())
            ):
                # The copies split symmetrically, so the midpoint cancels their O(sqrt(eps)) error.
                # Chained merges of 3+ close roots bias the location toward the last root.
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
    signs = _interior_signs(boundary, left, right)
    resolved_signs = signs[signs != 0]
    if resolved_signs.size == 0:
        _warn(
            f"sign of the boundary polynomial is unresolved throughout the segment "
            f"({left:.6g}, {right:.6g}) on the {boundary.label}; counting it as stable"
        )
        return True
    if np.any(resolved_signs != resolved_signs[0]):
        _warn(
            f"the boundary polynomial changes sign inside the segment ({left:.6g}, "
            f"{right:.6g}) on the {boundary.label}, so a crossing was missed; the "
            f"verdict for this segment is taken from its first resolved sample"
        )
    return bool(resolved_signs[0] < 0)


def _check_end_segment(signs: np.ndarray, stable: bool, message: str) -> None:
    """Warn if a resolved sample sign contradicts the verdict for an end segment."""
    bad = signs > 0 if stable else signs < 0
    if bad.any():
        _warn(message)


def _classify_stable_intervals(boundary: _Boundary) -> _AxisIntervals:
    """Group the segments between crossings into maximal stable runs, in the parameter t.

    Every verdict is based on the sign of B. The first segment is judged at the origin
    (:attr:`_Boundary.origin_stable`), the last by :attr:`_Boundary.far_stable`, and the rest by
    :func:`_segment_stable`. The first and last segments are also sampled, and the no-crossing case
    is checked, to warn about a missed crossing without changing the verdict.
    """
    crossings = _crossings(boundary)
    if not crossings:
        if boundary.origin_stable != boundary.far_stable:
            _warn(
                f"no crossings were found on the {boundary.label}, but the signs at the "
                f"origin and far from it disagree, so a crossing was missed (or lies at or "
                f"below the root-zero tolerance); the whole axis is "
                f"classified from the sign at the origin"
            )
        return _AxisIntervals([], boundary.origin_stable)  # uniformly stable or unstable

    _check_end_segment(
        _interior_signs(boundary, 0.0, crossings[0]),
        boundary.origin_stable,
        f"a resolved sample inside the first segment (0, {crossings[0]:.6g}) on the "
        f"{boundary.label} disagrees with the sign at the origin, so a crossing was missed "
        f"(or lies at or below the root-zero tolerance); the verdict for this segment is still "
        f"taken from the origin",
    )
    _check_end_segment(
        boundary.signs(crossings[-1] * _FAR_SAMPLE_FACTORS),
        boundary.far_stable,
        f"a resolved sample beyond the last crossing ({crossings[-1]:.6g}) on the "
        f"{boundary.label} disagrees with the sign far from the origin, so a crossing was "
        f"missed; the verdict for this segment is still taken from the leading coefficients",
    )

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


# ---------------------------------------------------------------------------
# Stiff properties
# ---------------------------------------------------------------------------


def _effective_degree(poly: np.poly1d, error: np.ndarray) -> int:
    """Degree of *poly* counting only coefficients resolved against their round-off bounds."""
    significant = np.flatnonzero(resolved(poly.c, error))
    return len(poly.c) - 1 - int(significant[0]) if significant.size else -1


def r_infinity(phi: StabilityFunction) -> float:
    """|R(infinity)| = lim |p(z) / q(z)| as |z| -> infinity.

    0.0 if deg p < deg q, |p_lead / q_lead| if the degrees are equal, and math.inf if deg p >
    deg q (every explicit method). Degrees count only resolved coefficients, so a leading
    coefficient that cancelled algebraically (TRBDF2, Billington) gives exactly 0.0.
    """
    p_deg = _effective_degree(phi.p, phi.p_error)
    q_deg = _effective_degree(phi.q, phi.q_error)
    if p_deg < q_deg:
        return 0.0
    if p_deg > q_deg:
        return math.inf
    return float(abs(phi.p.c[len(phi.p.c) - 1 - p_deg] / phi.q.c[len(phi.q.c) - 1 - q_deg]))


def _poles(phi: StabilityFunction) -> np.ndarray:
    """Non-removable poles of phi: roots of q at which p does not vanish (within round-off)."""
    q_deg = _effective_degree(phi.q, phi.q_error)
    if q_deg < 1:
        return np.array([], dtype=complex)
    poles = np.roots(phi.q.c[len(phi.q.c) - 1 - q_deg :])
    # A pole cancelled by a zero of p is removable: |p| there is unresolved from zero.
    bound = phi.p_error + phi.p.c.size * EPS * np.abs(phi.p.c)
    keep = [z for z in poles if abs(np.polyval(phi.p.c, z)) > 8.0 * np.polyval(bound, abs(z))]
    return np.array(keep, dtype=complex)


def _left_half_plane_poles(poles: np.ndarray) -> np.ndarray:
    """Poles with Re < -TOL.pole_rtol |pole|."""
    return poles[poles.real < -TOL.pole_rtol * np.abs(poles)]


def _a_stable(phi: StabilityFunction, imag_extent: float | None, poles: np.ndarray) -> bool:
    """A-stability from a precomputed imaginary-axis extent and the poles."""
    if imag_extent is not None:
        return False
    if _left_half_plane_poles(poles).size:
        return False
    if poles.size and np.any(np.abs(poles.real) <= TOL.pole_rtol * np.abs(poles)):
        _warn(
            "a pole of the stability function lies within the pole tolerance of the imaginary "
            "axis; A-stability is decided by the imaginary-axis classification alone"
        )
    return True


def is_a_stable(phi: StabilityFunction) -> bool:
    """True if |phi(z)| <= 1 for every z with Re z <= 0.

    Requires the imaginary axis to lie entirely in the region and no non-removable pole in the open
    left half-plane; by the maximum modulus principle that is sufficient. Poles within
    ``TOL.pole_rtol`` of the imaginary axis warn, and the imaginary-axis verdict decides.
    """
    return _a_stable(phi, axis_extent(phi, "imag"), _poles(phi))


def is_l_stable(phi: StabilityFunction) -> bool:
    """True if phi is A-stable and |phi(infinity)| = 0."""
    return is_a_stable(phi) and r_infinity(phi) == 0.0


def _a_alpha_angle(
    phi: StabilityFunction, real_extent: float | None, a_stable: bool
) -> float | None:
    """A(alpha) from precomputed real-axis extent and A-stability; see :func:`a_alpha_angle`."""
    if real_extent is not None:
        return None
    if a_stable:
        return 90.0

    # Upper bound: a pole in the open left half-plane stops the sector at its own angle.
    poles = _poles(phi)
    left = poles[poles.real < 0.0]
    cap = 90.0
    if left.size:
        cap = min(cap, float(np.degrees(np.arctan2(np.abs(left.imag), -left.real)).min()))

    # Guard scan at fixed steps below the cap; the first unstable angle ends the search window.
    lo, hi = 0.0, cap
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", StabilityWarning)
        step = _A_ALPHA_SCAN_STEP_DEG
        for angle in np.arange(step, cap, step):
            if _ray_stable(phi, float(angle)):
                lo = float(angle)
            else:
                hi = float(angle)
                break
        while hi - lo > TOL.angle:
            mid = 0.5 * (lo + hi)
            if _ray_stable(phi, mid):
                lo = mid
            else:
                hi = mid
    return lo


def a_alpha_angle(phi: StabilityFunction) -> float | None:
    """The A(alpha) angle in degrees: the largest alpha with the sector |arg(-z)| <= alpha stable.

    Returns None if the negative real axis is not entirely in the region (every explicit method,
    and any |R(infinity)| > 1), 90.0 if and only if phi is A-stable, and 0.0 if the real axis is
    stable but no sector around it is. The search is exact per ray (the same roots, tolerances and
    closed-region convention as the axes) and then bisects the angle to ``TOL.angle``. It assumes a
    consistent method, phi(z) = 1 + z + O(z^2); see the module notes for why the stable angles form
    one interval.
    """
    real_extent = axis_extent(phi, "real")
    if real_extent is not None:
        return None
    return _a_alpha_angle(phi, real_extent, is_a_stable(phi))


def threshold_factor(phi: StabilityFunction) -> float | None:
    """Linear threshold factor of a polynomial phi; None if phi is rational.

    The largest r such that phi(z) = sum_k g_k (1 + z/r)^k with every g_k >= 0, equivalently every
    Taylor coefficient of phi at z = -r is nonnegative. It is at least the SSP coefficient. It is
    0.0 if some coefficient of p is not positive. Implicit methods return None (not computed).
    """
    if _effective_degree(phi.q, phi.q_error) > 0:
        return None
    degree = _effective_degree(phi.p, phi.p_error)
    if degree < 1:
        return None
    coeffs = np.asarray(phi.p.c[len(phi.p.c) - 1 - degree :], float)[::-1]  # ascending
    errors = np.asarray(phi.p_error[len(phi.p.c) - 1 - degree :], float)[::-1]
    if not np.all(resolved(coeffs, errors) & (coeffs > 0.0)):
        return 0.0

    ks = np.arange(degree + 1)
    binom = np.array([[math.comb(j, k) if j >= k else 0 for j in ks] for k in ks], float)

    def feasible(r: float) -> bool:
        # g_k(r) = sum_{j >= k} p_j C(j, k) (-r)^(j - k)
        powers = np.where(
            ks[None, :] >= ks[:, None], r ** np.maximum(ks[None, :] - ks[:, None], 0), 0.0
        )
        terms = binom * coeffs[None, :] * powers
        signed = terms * (-1.0) ** np.maximum(ks[None, :] - ks[:, None], 0)
        slack = 16.0 * EPS * terms.sum(axis=1)
        return bool(np.all(signed.sum(axis=1) >= -slack))

    # g_{n-1}(r) = p_{n-1} - n p_n r >= 0 bounds r above; the feasible set is [0, r*].
    lo, hi = 0.0, float(coeffs[degree - 1] / (degree * coeffs[degree]))
    if feasible(hi):
        return hi
    while hi - lo > 4.0 * EPS * hi:
        mid = 0.5 * (lo + hi)
        if feasible(mid):
            lo = mid
        else:
            hi = mid
    return lo


def order_star_function(phi: StabilityFunction, Z) -> np.ndarray:
    """F(z) = log|p(z)| - log|q(z)| - Re z, the log of |phi(z)| / |e^z|, for scalar or array Z.

    The order star is {F > 0}. The exponential is folded into the subtraction so it never
    overflows. NaN occurs only at exact common zeros of p and q; inf at a pole, -inf at a zero.
    """
    Z = np.asarray(Z, dtype=complex)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log(np.abs(phi.p(Z))) - np.log(np.abs(phi.q(Z))) - Z.real


__all__ = [
    "TOL",
    "StabilityFunction",
    "StabilityWarning",
    "a_alpha_angle",
    "axis_extent",
    "imag_stable_intervals",
    "is_a_stable",
    "is_l_stable",
    "max_axis_crossing",
    "order_star_function",
    "r_infinity",
    "stability_magnitude",
    "threshold_factor",
]
