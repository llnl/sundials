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
# Module to find and setup PETSC.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

include_guard(GLOBAL)

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# Using PETSc requires building with MPI enabled
if(SUNDIALS_ENABLE_PETSC AND NOT SUNDIALS_ENABLE_MPI)
  message(
    FATAL_ERROR
      "MPI is required for PETSc support. Set SUNDIALS_ENABLE_MPI to ON.")
endif()

if(SUNDIALS_PRECISION MATCHES "EXTENDED")
  message(
    FATAL_ERROR
      "SUNDIALS is not compatible with PETSc when using ${SUNDIALS_PRECISION} precision"
  )
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

find_package(PETSC 3.5.0 REQUIRED)

message(STATUS "PETSC_DIR:          ${PETSC_DIR}")
message(STATUS "PETSC_ARCH:         ${PETSC_ARCH}")
message(STATUS "PETSC_LIBRARIES:    ${PETSC_LIBRARIES}")
message(STATUS "PETSC_INCLUDE_DIRS: ${PETSC_INCLUDE_DIRS}")
message(STATUS "PETSC_INDEX_SIZE:   ${PETSC_INDEX_SIZE}")
message(STATUS "PETSC_PRECISION:    ${PETSC_PRECISION}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

if(SUNDIALS_ENABLE_PETSC_CHECKS)
  # No need for any compile tests because the FindPETSC module does compile
  # tests already.

  message(CHECK_START "Testing PETSc")

  if(NOT ("${SUNDIALS_INDEX_SIZE}" MATCHES "${PETSC_INDEX_SIZE}"))
    message(CHECK_FAIL "failed")
    string(
      CONCAT _err_msg_string
             "PETSc not functional due to index size mismatch:\n"
             "SUNDIALS_INDEX_SIZE=${SUNDIALS_INDEX_SIZE}, "
             "but PETSc was built with ${PETSC_INDEX_SIZE}-bit indices\n"
             "PETSC_DIR: ${PETSC_DIR}\n")
    message(FATAL_ERROR "${_err_msg_string}")
  endif()

  string(TOUPPER "${PETSC_PRECISION}" _petsc_precision)
  string(TOUPPER "${SUNDIALS_PRECISION}" _sundials_precision)
  if(NOT ("${_sundials_precision}" MATCHES "${_petsc_precision}"))
    message(CHECK_FAIL "failed")
    string(
      CONCAT _err_msg_string
             "PETSc not functional due to real type precision mismatch:\n"
             "SUNDIALS_PRECISION=${_sundials_precision}, "
             "but PETSc was built with ${_petsc_precision} precision\n"
             "PETSC_DIR: ${PETSC_DIR}\n")
    message(FATAL_ERROR "${_err_msg_string}")
  endif()

  message(CHECK_PASS "success")

else()
  message(STATUS "Skipped PETSC checks.")
endif()
