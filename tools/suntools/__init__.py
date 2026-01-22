#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025, Lawrence Livermore National Security,
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

"""
This is a Python library of utilities SUNDIALS developer may find useful.
Right now it consists of the following modules:

- `logs`: this module has functions for parsing logs produced by `SUNLogger`.
- `table`: this module has functions for parsing stats output by SUNDIALS
           integrators in the table format.
- `csv`: this module has functions for parsing stats output by SUNDIALS
         integrators in the CSV format.
"""
