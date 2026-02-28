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
# Module to find and setup SuperLU_DIST correctly.
# Created from the SundialsTPL.cmake template.
# All SUNDIALS modules that find and setup a TPL must:
#
# 1. Check to make sure the SUNDIALS configuration and the TPL is compatible.
# 2. Find the TPL.
# 3. Check if the TPL works with SUNDIALS, UNLESS the override option
# <TPL>_WORKS is TRUE - in this case the tests should not be performed and it
# should be assumed that the TPL works with SUNDIALS.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

if(NOT DEFINED SUNDIALS_SUPERLUDIST_INCLUDED)
  set(SUNDIALS_SUPERLUDIST_INCLUDED)
else()
  return()
endif()

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
if(ENABLE_SUPERLUDIST AND NOT ENABLE_MPI)
  message(
    FATAL_ERROR
      "MPI is required for SuperLU DIST support. Set ENABLE_MPI to ON.")
endif()

# Using SUPERLUDIST with OpenMP requires building with OpenMP enabled
if(ENABLE_SUPERLUDIST
   AND SUPERLUDIST_OpenMP
   AND NOT ENABLE_OPENMP)
  message(
    FATAL_ERROR
      "OpenMP is required for SuperLU DIST support. Set ENABLE_OPENMP to ON.")
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
if(SUPERLUDIST_FOUND AND (NOT SUPERLUDIST_WORKS))

  if(SUPERLUDIST_CUDA AND (NOT ENABLE_CUDA))
    message(
      FATAL_ERROR
        "SuperLU_DIST was built with CUDA but SUNDIALS does not have CUDA enabled. Set ENABLE_CUDA=TRUE."
    )
  endif()

  if(SUPERLUDIST_HIP AND (NOT ENABLE_HIP))
    message(
      FATAL_ERROR
        "SuperLU_DIST was built with HIP but SUNDIALS does not have HIP enabled. Set ENABLE_HIP=TRUE."
    )
  endif()

  # Check index size
  if(NOT (SUNDIALS_INDEX_SIZE STREQUAL SUPERLUDIST_INDEX_SIZE))
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

  message(STATUS "Checking if SuperLU_DIST works with SUNDIALS... OK")
  set(SUPERLUDIST_WORKS
      TRUE
      CACHE BOOL "SuperLU_DIST works with SUNDIALS as configured" FORCE)

elseif(SUPERLUDIST_FOUND AND SUPERLUDIST_WORKS)
  message(
    STATUS
      "Skipped SuperLU_DIST tests, assuming SuperLU_DIST works with SUNDIALS. Set SUPERLUDIST_WORKS=FALSE to (re)run compatibility test."
  )
endif()
