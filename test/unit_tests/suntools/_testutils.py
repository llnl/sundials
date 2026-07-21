#!/usr/bin/env python3
# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos
# ---------------------------------------------------------------
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
# ---------------------------------------------------------------

import os
import sys


def add_repo_suntools_to_path():
    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
    suntools_src = os.path.join(repo_root, "suntools", "src")
    if suntools_src not in sys.path:
        sys.path.insert(0, suntools_src)
