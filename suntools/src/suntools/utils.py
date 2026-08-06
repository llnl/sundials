#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# Shared suntools utilities.
# -----------------------------------------------------------------------------

from __future__ import annotations

import numpy as np

# Multiplier on machine epsilon separating round-off from genuine values.
_ROUNDOFF_FACTOR = 256.0

# Safety factor on a propagated error bound.
_SIGNIFICANCE_FACTOR = 8.0

EPS = float(np.finfo(float).eps)


def str2num(s):
    """Try to convert a string to an int or float"""

    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s


def resolved(values, bounds) -> np.ndarray:
    """Which entries are distinguishable from their error bound; scalar or array."""
    return np.abs(values) > _SIGNIFICANCE_FACTOR * np.asarray(bounds, float)


def roundoff_tol(values, factor: float = _ROUNDOFF_FACTOR) -> float:
    """Magnitude below which an entry of *values* is indistinguishable from zero."""
    values = np.asarray(values)
    scale = max(1.0, float(np.max(np.abs(values)))) if values.size else 1.0
    return factor * EPS * scale
