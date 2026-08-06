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

"""Runge-Kutta Butcher tables and linear stability functions."""

from __future__ import annotations

from functools import cached_property

import numpy as np

from .rk_stability import StabilityFunction
from .utils import EPS, resolved, roundoff_tol


def _det_poly(build_matrix, stages: int):
    """Polynomial coefficients of det(M(z)) for a real-affine matrix M(z).

    Since M(z) is affine, det(M(z)) is a polynomial in z of degree at most ``stages``, and because
    A and b are real the polynomial coefficients are real. Writing

        f(z) = sum_j a_j z^j for j = 0,...,stages

    and evaluating at any n roots of unity, z_k = e^{2 pi i k / n}, gives

        f(z_k) = sum_j a_j e^{2 pi i j k / n}.

    Multiplying by e^{-2 pi i m k / n}, summing over k, and using orthogonality of the roots of
    unity gives a_m = (1/n) FFT(f)[m], where a_m = 0 for m > stages. We use twice as many values as
    necessary, n = 2(stages + 1), to estimate the roundoff error. Coefficients are returned in
    descending order with negligible leading coefficients dropped.
    """
    n_nodes = 2 * (stages + 1)
    nodes = np.exp(2j * np.pi * np.arange(n_nodes) / n_nodes)
    det_values = np.array([np.linalg.det(build_matrix(z)) for z in nodes])
    coeffs = np.fft.fft(det_values) / n_nodes  # coeffs[j] = coefficient of z^j
    # The error in computing the determinant from an LU factorization with partial pivoting scales
    # with ||M||_2 ||adj(M)||_2 = sv[0] * prod(sv[:-1]). This bound assumes worst-case pivot growth
    # and worst-case alignment of the backward error which can be enough to trim real leading
    # coefficients. Instead we over sample and use coefficients that should vanish exactly (index
    # above `stages`), to measure the noise floor.
    safety = 8.0
    noise = safety * float(np.max(np.abs(coeffs[stages + 1:])))
    error = np.full(stages + 1, max(noise, EPS * float(np.max(np.abs(det_values)))))
    if np.any(resolved(coeffs[:stages + 1].imag, error)):
        raise ValueError("det(M(z)) recovered a non-negligible imaginary part.")
    return _trim_with_error(coeffs[stages::-1].real, error)


def _trim_with_error(coefficients, error):
    """Drop leading coefficients that are round-off compared to the error level."""
    coefficients = np.asarray(coefficients, float)
    error = np.asarray(error, float)
    significant = np.nonzero(resolved(coefficients, error))[0]
    if significant.size == 0:
        # Cannot happen for a consistent method, whose p and q both have constant term 1.
        return np.zeros(1), error[-1:]
    top = int(significant[0])
    return coefficients[top:], error[top:]


def _pq_from_weights(A, weights, stages, explicit):
    """Build p(z) and q(z) in phi(z) = p(z) / q(z).

        phi(z) = det(I - zA + z e b^T) / det(I - z A) = 1 + z b^T (I - z A)^{-1} e

    where A holds the stage coefficients, b the solution or embedding weights, and e is a vector of
    ones. For explicit methods q(z) = 1 and phi is a polynomial whose coefficients follow from a
    recurrence relation; for implicit methods determinants are evaluated with :func:`_det_poly`.
    """
    if explicit:
        coeffs = np.zeros(stages + 1)
        error = np.zeros(stages + 1)
        abs_A, abs_b = np.abs(A), np.abs(weights)
        coeffs[0] = 1.0
        # A is strictly lower triangular, so det(I - zA) = 1 and (I - zA)^{-1} = sum_k z^k A^k
        # terminates: a_0 = 1, a_{k+1} = b^T w_k with w_k = A^k e built as w_{k+1} = A w_k.
        Ak_ones = np.ones(stages)
        # W_k = |A|^k e bounds |w_k| and its error E_k, via E_{k+1} = |A| E_k for existing error
        # plus (s + 1) EPS |A| W_k for fresh rounding and the tableau's own EPS storage error.
        magnitude = np.ones(stages)
        bound = np.zeros(stages)
        for k in range(stages):
            coeffs[k + 1] = weights @ Ak_ones
            error[k + 1] = abs_b @ bound + (stages + 1) * EPS * (abs_b @ magnitude)
            Ak_ones = A @ Ak_ones
            bound = abs_A @ bound + (stages + 1) * EPS * (abs_A @ magnitude)
            magnitude = abs_A @ magnitude
        p, p_error = _trim_with_error(coeffs[::-1], error[::-1])
        return np.poly1d(p), np.poly1d([1.0]), p_error, np.zeros(1)

    I = np.eye(stages)
    ones = np.ones(stages)
    p, p_error = _det_poly(lambda z: I - z * A + z * np.outer(ones, weights), stages)
    q, q_error = _det_poly(lambda z: I - z * A, stages)
    return np.poly1d(p), np.poly1d(q), p_error, q_error


class ButcherTable:
    """Runge-Kutta Butcher table.

    Besides the tableau arrays (A, b, c) this class builds the method's linear stability function,
    see :class:`~suntools.rk_stability.StabilityFunction`.

    Attributes
    ----------
    name : str                               method name (e.g. "ARKODE_DORMAND_PRINCE_7_4_5")
    stages : int                             number of stages
    A : (stages, stages) ndarray             stage coefficient matrix
    b : (stages,) ndarray                    solution weights
    c : (stages,) ndarray                    abscissae
    b_embedded : (stages,) ndarray or None   embedding weights
    method_order : int or None               order of the main method
    embedding_order : int or None            order of the embedded method
    """

    def __init__(
        self, name, stages, A, b, c, b_embedded=None, method_order=None, embedding_order=None
    ):
        self.name = name
        self.stages = stages
        self.A = A
        self.b = b
        self.c = c
        self.b_embedded = b_embedded
        self.method_order = method_order
        self.embedding_order = embedding_order

    @cached_property
    def is_explicit(self) -> bool:
        """True if A is strictly lower triangular"""
        return not bool(np.any(np.abs(np.triu(self.A)) > roundoff_tol(self.A)))

    @property
    def has_embedding(self) -> bool:
        """True if there is an embedded method."""
        return self.b_embedded is not None

    def kind(self) -> str:
        """Classify the method from the structure of A."""
        if self.is_explicit:
            return "explicit (ERK)"
        tolerance = roundoff_tol(self.A)
        if np.any(np.abs(np.triu(self.A, 1)) > tolerance):
            return "fully implicit (IRK)"
        diagonal = np.diag(self.A)
        if abs(diagonal[0]) <= tolerance:
            rest = diagonal[1:]
            if np.all(np.abs(rest - rest[0]) <= tolerance):
                return "explicit first stage, singly diagonally implicit (ESDIRK)"
            return "explicit first stage, diagonally implicit (EDIRK)"
        if np.all(np.abs(diagonal - diagonal[0]) <= tolerance):
            return "singly diagonally implicit (SDIRK)"
        return "diagonally implicit (DIRK)"

    @cached_property
    def stability_function(self) -> StabilityFunction:
        """The method's linear stability function."""
        return StabilityFunction(
            *_pq_from_weights(self.A, self.b, self.stages, self.is_explicit),
            order=self.method_order,
        )

    @cached_property
    def embedded_stability_function(self) -> StabilityFunction | None:
        """The embedding's linear stability function, or None if no embedding."""
        if self.b_embedded is None:
            return None
        return StabilityFunction(
            *_pq_from_weights(self.A, self.b_embedded, self.stages, self.is_explicit),
            order=self.embedding_order,
        )

    def __str__(self) -> str:
        """Format the tableau for print(table)."""
        fmt = "{:.6g}".format  # Table entries are printed with this format
        rows = [[fmt(v) for v in row] for row in self.A]
        crow = [fmt(v) for v in self.c]
        brow = [fmt(v) for v in self.b]
        drow = [fmt(v) for v in self.b_embedded] if self.has_embedding else None

        width = max(len(v) for v in [v for row in rows for v in row] + brow + (drow or ["0"]))
        cwidth = max(len(v) for v in crow)

        def body(values):
            return " ".join(v.rjust(width) for v in values)

        metadata = []
        if self.method_order is not None:
            metadata.append(f"order q = {self.method_order}")
        if self.embedding_order is not None:
            metadata.append(f"embedding order p = {self.embedding_order}")
        metadata += [f"stages s = {self.stages}", self.kind()]

        lines = [f"Butcher table: {self.name}", "=" * (len(self.name) + 15), "  ".join(metadata), ""]
        lines += [f"{crow[i].rjust(cwidth)} | {body(rows[i])}" for i in range(self.stages)]
        lines.append("-" * (cwidth + 3 + self.stages * (width + 1)))
        lines.append(f"{'':>{cwidth}} | {body(brow)}   (b)")
        if drow is not None:
            lines.append(f"{'':>{cwidth}} | {body(drow)}   (b-hat, embedding)")
        return "\n".join(lines)

    def __repr__(self):
        """Compact one-line summary."""
        parts = [
            f"name={self.name!r}",
            f"stages={self.stages}",
            f"method_order={self.method_order}",
        ]
        if self.has_embedding:
            parts.append(f"embedding_order={self.embedding_order}")
        return f"ButcherTable({', '.join(parts)})"
