# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------
# nonlinear solver level CMakeLists.txt for SUNDIALS
# ------------------------------------------------------------------------------

# required modules
add_subdirectory(newton)
add_subdirectory(fixedpoint)

if(SUNDIALS_ENABLE_SUNNONLINSOL_PETSCSNES)
  add_subdirectory(petscsnes)
endif()
