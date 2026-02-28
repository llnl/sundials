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
# Module to find and setup MAGMA correctly.
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

if(NOT DEFINED SUNDIALS_MAGMA_INCLUDED)
  set(SUNDIALS_MAGMA_INCLUDED)
else()
  return()
endif()

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

if(SUNDIALS_PRECISION MATCHES "extended")
  message(
    FATAL_ERROR
      "SUNDIALS MAGMA interface is not compatible with extended precision")
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

find_package(MAGMA REQUIRED)

message(STATUS "MAGMA_VERSION:           ${MAGMA_VERSION}")
message(STATUS "MAGMA_LIBRARIES:         ${MAGMA_LIBRARIES}")
message(STATUS "MAGMA_INCLUDE_DIR:       ${MAGMA_INCLUDE_DIR}")
message(STATUS "SUNDIALS_MAGMA_BACKENDS: ${SUNDIALS_MAGMA_BACKENDS}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

if(MAGMA_FOUND AND (NOT MAGMA_WORKS))
  if(SUNDIALS_MAGMA_BACKENDS MATCHES "CUDA" AND NOT ENABLE_CUDA)
    message(
      FATAL_ERROR
        "SUNDIALS_MAGMA_BACKENDS includes CUDA but CUDA is not enabled. Set ENABLE_CUDA=ON or change the backend."
    )
  endif()
  if(SUNDIALS_MAGMA_BACKENDS MATCHES "HIP" AND NOT ENABLE_HIP)
    message(
      FATAL_ERROR
        "SUNDIALS_MAGMA_BACKENDS includes HIP but HIP is not enabled. Set ENABLE_HIP=ON or change the backend."
    )
  endif()

  set(MAGMA_WORKS
      TRUE
      CACHE BOOL "MAGMA works with SUNDIALS as configured" FORCE)
elseif(MAGMA_FOUND AND MAGMA_WORKS)
  message(STATUS "Skipped MAGMA tests, assuming MAGMA works with SUNDIALS.")
endif()
