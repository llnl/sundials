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
# test/unit_tests/sunnonlinsol level CMakeLists.txt for SUNDIALS
# ------------------------------------------------------------------------------

if(SUNDIALS_ENABLE_CALIPER AND SUNDIALS_BUILD_WITH_PROFILING)
  include_directories(${caliper_INCLUDE_DIR})
  set(EXE_EXTRA_LINK_LIBS ${EXE_EXTRA_LINK_LIBS} caliper)
endif()

# Always add the Newton and fixed point examples
add_subdirectory(newton)
add_subdirectory(fixedpoint)

if(BUILD_SUNNONLINSOL_PETSCSNES)
  add_subdirectory(petsc)
endif()
