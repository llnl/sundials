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

Given a method's linear stability function phi(z) = p(z)/q(z), this module describes the absolute
stability region { z : |phi(z)| <= 1 }.

Public API
----------
``StabilityFunction``            phi(Z) = p(Z)/q(Z).
``stability_magnitude(phi, Z)``  |phi(Z)|.

A plot draws the region as the |phi| = 1 contour of a grid evaluation, so |phi| is all it needs
from this module; nothing here solves for where the boundary sits.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class StabilityFunction:
    """phi(z) = p(z) / q(z)"""

    p: np.poly1d
    q: np.poly1d
    order: int | None = None  # RK method order


def stability_magnitude(phi, Z):
    """|phi(Z)| = |p(Z) / q(Z)| for scalar or array Z."""
    # Divide and invalid warnings are suppressed to keep batch plotting quiet when poles are
    # present; the values (including inf/nan) are preserved for plotting contours.
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.abs(phi.p(Z) / phi.q(Z))


__all__ = ["StabilityFunction", "stability_magnitude"]
