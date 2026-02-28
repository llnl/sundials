# ------------------------------------------------------------------------------
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
# Module to find and setup Kokkos correctly.
# Created from the SundialsTPL.cmake template.
# All SUNDIALS modules that find and setup a TPL must:
#
# 1. Check to make sure the SUNDIALS configuration and the TPL is compatible.
# 2. Find the TPL.
# 3. Check if the TPL works with SUNDIALS, UNLESS the override option
# TPL_WORKS is TRUE - in this case the tests should not be performed and it
# should be assumed that the TPL works with SUNDIALS.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

if(NOT DEFINED SUNDIALS_KOKKOS_INCLUDED)
  set(SUNDIALS_KOKKOS_INCLUDED)
else()
  return()
endif()

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------
find_package(Kokkos REQUIRED HINTS "${Kokkos_DIR}" NO_DEFAULT_PATH)

# We should be able to use Kokkos_DEVICES directly but it seems to get removed
# or unset in some CMake versions
set(KOKKOS_EXAMPLES_BACKENDS
    "${Kokkos_DEVICES}"
    CACHE STRING "Kokkos backends to build examples with")
mark_as_advanced(FORCE KOKKOS_EXAMPLES_BACKENDS)
message(STATUS "Kokkos VERSION: ${Kokkos_VERSION}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

if(Kokkos_FOUND AND (NOT KOKKOS_WORKS))
  message(STATUS "Checking if Kokkos works... OK")
  set(KOKKOS_WORKS
      TRUE
      CACHE BOOL "Kokkos works with SUNDIALS as configured" FORCE)
elseif(Kokkos_FOUND AND KOKKOS_WORKS)
  message(STATUS "Skipped Kokkos tests, assuming Kokkos works with SUNDIALS.")
endif()
