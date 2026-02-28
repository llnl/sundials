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
find_package(KokkosKernels REQUIRED HINTS "${KokkosKernels_DIR}"
             NO_DEFAULT_PATH)

message(STATUS "Kokkos Kernels VERSION: ${KokkosKernels_VERSION}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

if(KokkosKernels_FOUND AND (NOT KOKKOS_KERNELS_WORKS))
  message(STATUS "Checking if Kokkos Kernels works... OK")
  set(KOKKOS_KERNELS_WORKS
      TRUE
      CACHE BOOL "Kokkos Kernels works with SUNDIALS as configured" FORCE)
elseif(KokkosKernels_FOUND AND KOKKOS_WORKS)
  message(
    STATUS
      "Skipped Kokkos Kernels tests, assuming Kokkos Kernels works with SUNDIALS."
  )
endif()
