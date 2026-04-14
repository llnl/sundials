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
# Module to find and setup Trilinos.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

include_guard(GLOBAL)

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

# Find Trilinos
find_package(
  Trilinos REQUIRED
  COMPONENTS Tpetra HINTS "${Trilinos_DIR}/lib/cmake/Trilinos"
             "${Trilinos_DIR}")

message(STATUS "Trilinos Libraries: ${Trilinos_LIBRARIES}")
message(STATUS "Trilinos Includes: ${Trilinos_INCLUDE_DIRS}")
message(STATUS "Trilinos Devices: ${Kokkos_DEVICES}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

# Does not currently work with Trilinos imported targets due to an error from
# evaluating generator expression: $<LINK_LANGUAGE:CXX> may only be used with
# binary targets to specify link libraries, link directories, link options and
# link depends.
