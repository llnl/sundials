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
# Module to find and setup SUPERLUMT.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

include_guard(GLOBAL)

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# SUPERLUMT does not support extended precision
if(SUNDIALS_PRECISION MATCHES "EXTENDED")
  message(
    FATAL_ERROR
      "SUPERLUMT is not compatible with ${SUNDIALS_PRECISION} precision")
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

find_package(SUPERLUMT REQUIRED)

message(STATUS "SUPERLUMT_LIBRARIES:   ${SUPERLUMT_LIBRARIES}")
message(STATUS "SUPERLUMT_INCLUDE_DIR: ${SUPERLUMT_INCLUDE_DIR}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

if(SUNDIALS_ENABLE_SUPERLUMT_CHECKS)

  message(CHECK_START "Testing SuperLU_MT")

  # Create the test directory
  set(TEST_DIR ${PROJECT_BINARY_DIR}/SUPERLUMT_TEST)

  # Create a C source file which uses SUPERLUMT types
  file(
    WRITE ${TEST_DIR}/test.c
    "\#include \"slu_mt_ddefs.h\"\n"
    "int main(void) {\n"
    "SuperMatrix *A;\n"
    "NCformat *Astore;\n"
    "A = NULL;\n"
    "Astore = NULL;\n"
    "if (A != NULL || Astore != NULL) return(1);\n"
    "else return(0);\n"
    "}\n")

  # Attempt to build and link the test executable, pass --debug-trycompile to
  # the cmake command to save build files for debugging
  try_compile(
    COMPILE_OK ${TEST_DIR}
    ${TEST_DIR}/test.c
    LINK_LIBRARIES SUNDIALS::SUPERLUMT
    OUTPUT_VARIABLE COMPILE_OUTPUT)

  # Check the result
  if(COMPILE_OK)
    message(CHECK_PASS "success")
  else()
    message(CHECK_FAIL "failed")
    file(WRITE ${TEST_DIR}/compile.out "${COMPILE_OUTPUT}")
    message(
      FATAL_ERROR
        "Could not compile SuperLU_MT test. Check output in ${TEST_DIR}/compile.out"
    )
  endif()

else()
  message(STATUS "Skipped SuperLU_MT checks.")
endif()
