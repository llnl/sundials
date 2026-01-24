# -----------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2025, Lawrence Livermore National Security,
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
# -----------------------------------------------------------------
# Base classes for problems
# -----------------------------------------------------------------


class ODEProblem:

    def set_init_cond(self, y0vec):
        raise NotImplementedError("Subclasses must implement the set_init_cond method.")

    def solution(self, y0vec, yvec, t):
        raise NotImplementedError("Subclasses must implement the set_init_cond method.")
