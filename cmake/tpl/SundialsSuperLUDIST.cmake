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
# Module to find and setup SuperLU_DIST.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

include_guard(GLOBAL)

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# SuperLU_DIST only supports double precision
if(SUNDIALS_PRECISION MATCHES "SINGLE" OR SUNDIALS_PRECISION MATCHES "EXTENDED")
  message(
    FATAL_ERROR
      "SuperLU_DIST is not compatible with ${SUNDIALS_PRECISION} precision")
endif()

# Using SUPERLUDIST requires building with MPI enabled
if(SUNDIALS_ENABLE_SUPERLUDIST AND NOT SUNDIALS_ENABLE_MPI)
  message(
    FATAL_ERROR
      "MPI is required for SuperLU DIST support. Set SUNDIALS_ENABLE_MPI to ON."
  )
endif()

# Using SUPERLUDIST with OpenMP requires building with OpenMP enabled
if(SUNDIALS_ENABLE_SUPERLUDIST
   AND SUPERLUDIST_OpenMP
   AND NOT SUNDIALS_ENABLE_OPENMP)
  message(
    FATAL_ERROR
      "OpenMP is required for SuperLU DIST support. Set SUNDIALS_ENABLE_OPENMP to ON."
  )
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

# We need MPI for SuperLU_DIST support
include(SundialsMPI)

# Try to find SuperLU_DIST
find_package(SUPERLUDIST 7.0.0 REQUIRED)

message(STATUS "SUPERLUDIST_VERSION:        ${SUPERLUDIST_VERSION}")
message(STATUS "SUPERLUDIST_LINK_LIBRARIES: ${SUPERLUDIST_LINK_LIBRARIES}")
message(STATUS "SUPERLUDIST_INCLUDE_DIRS:   ${SUPERLUDIST_INCLUDE_DIRS}")
message(STATUS "SUPERLUDIST_INDEX_SIZE:     ${SUPERLUDIST_INDEX_SIZE}")
message(STATUS "SUPERLUDIST_OpenMP:         ${SUPERLUDIST_OpenMP}")
message(STATUS "SUPERLUDIST_CUDA:           ${SUPERLUDIST_CUDA}")
message(STATUS "SUPERLUDIST_ROCM:           ${SUPERLUDIST_ROCM}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

# If we have the SuperLU_DIST libraries, test them
if(SUNDIALS_ENABLE_SUPERLUDIST_CHECKS)

  message(CHECK_START "Testing SuperLU_DIST")

  if(SUPERLUDIST_CUDA AND (NOT SUNDIALS_ENABLE_CUDA))
    message(CHECK_FAIL "failed")
    message(
      FATAL_ERROR
        "SuperLU_DIST was built with CUDA but SUNDIALS does not have CUDA enabled. Set SUNDIALS_ENABLE_CUDA=TRUE."
    )
  endif()

  if(SUPERLUDIST_ROCM AND (NOT SUNDIALS_ENABLE_HIP))
    message(CHECK_FAIL "failed")
    message(
      FATAL_ERROR
        "SuperLU_DIST was built with HIP but SUNDIALS does not have HIP enabled. Set SUNDIALS_ENABLE_HIP=TRUE."
    )
  endif()

  # Check index size
  if(NOT (SUNDIALS_INDEX_SIZE STREQUAL SUPERLUDIST_INDEX_SIZE))
    message(CHECK_FAIL "failed")
    set(_err_msg_string
        "SuperLU_DIST not functional due to index size mismatch:\n")
    string(
      APPEND
      _err_msg_string
      "SUNDIALS_INDEX_SIZE=${SUNDIALS_INDEX_SIZE}, but SuperLU_DIST was built with ${SUPERLUDIST_INDEX_SIZE}-bit indices\n"
    )
    string(APPEND _err_msg_string
           "SUPERLUDIST_INCLUDE_DIRS: ${SUPERLUDIST_INCLUDE_DIRS}\n")
    message(FATAL_ERROR "${_err_msg_string}")
  endif()

  message(CHECK_PASS "success")

else()
  message(STATUS "Skipped SuperLU_DIST checks.")
endif()
