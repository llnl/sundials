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
# Module to find and setup MAGMA.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

include_guard(GLOBAL)

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

if(SUNDIALS_ENABLE_MAGMA_CHECKS)

  message(CHECK_START "Testing MAGMA")

  if(SUNDIALS_MAGMA_BACKENDS MATCHES "CUDA" AND NOT SUNDIALS_ENABLE_CUDA)
    message(CHECK_FAIL "failed")
    message(
      FATAL_ERROR
        "SUNDIALS_MAGMA_BACKENDS includes CUDA but CUDA is not enabled. Set SUNDIALS_ENABLE_CUDA=ON or change the backend."
    )
  endif()

  if(SUNDIALS_MAGMA_BACKENDS MATCHES "HIP" AND NOT SUNDIALS_ENABLE_HIP)
    message(CHECK_FAIL "failed")
    message(
      FATAL_ERROR
        "SUNDIALS_MAGMA_BACKENDS includes HIP but HIP is not enabled. Set SUNDIALS_ENABLE_HIP=ON or change the backend."
    )
  endif()

  message(CHECK_PASS "success")

else()
  message(STATUS "Skipped MAGMA checks.")
endif()
