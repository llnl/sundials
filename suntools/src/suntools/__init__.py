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

"""Utilities for working with SUNDIALS output and applications.

The package contains the following modules:

- ``logs``
   Functions for parsing and filtering logs produced by :c:type:`SUNLogger`.
- ``table``
   Functions for parsing statistics written in SUNDIALS table format.
- ``csv``
   Functions for parsing statistics written in SUNDIALS CSV format.
- ``tune``
   Configuration models and runners for tuning SUNDIALS applications by
   appending ``SetOptions`` parameters to an executable command.
- `rk`: this namespace contains Runge-Kutta table parsing, stability, and
       plotting utilities.
"""
