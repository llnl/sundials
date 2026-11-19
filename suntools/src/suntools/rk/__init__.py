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
from .plotting import PlotOptions, plot_stability_region
from .stability import StabilityFunction, stability_magnitude
from .table_parser import parse_butcher_tables

__all__ = [
    "ButcherTable",
    "PlotOptions",
    "StabilityFunction",
    "parse_butcher_tables",
    "plot_stability_region",
    "stability_magnitude",
]
