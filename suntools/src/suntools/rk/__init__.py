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

"""Runge-Kutta utilities for suntools."""

from .butcher_table import ButcherTable
from .metrics import (
    StabilityMetrics,
    format_metrics_csv,
    format_metrics_rst,
    ssp_coefficient,
    stability_metrics,
)
from .plotting import PlotOptions, plot_order_star, plot_stability_region
from .stability import (
    TOL,
    StabilityFunction,
    StabilityWarning,
    a_alpha_angle,
    axis_extent,
    imag_stable_intervals,
    is_a_stable,
    is_l_stable,
    max_axis_crossing,
    order_star_function,
    r_infinity,
    stability_magnitude,
    threshold_factor,
)
from .table_parser import parse_butcher_tables

__all__ = [
    "TOL",
    "ButcherTable",
    "PlotOptions",
    "StabilityFunction",
    "StabilityMetrics",
    "StabilityWarning",
    "a_alpha_angle",
    "axis_extent",
    "format_metrics_csv",
    "format_metrics_rst",
    "imag_stable_intervals",
    "is_a_stable",
    "is_l_stable",
    "max_axis_crossing",
    "order_star_function",
    "parse_butcher_tables",
    "plot_order_star",
    "plot_stability_region",
    "r_infinity",
    "ssp_coefficient",
    "stability_magnitude",
    "stability_metrics",
    "threshold_factor",
]
