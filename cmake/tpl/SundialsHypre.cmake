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
# Module to find and setup HYPRE.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

include_guard(GLOBAL)

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# Using hypre requires building with MPI enabled
if(NOT SUNDIALS_ENABLE_MPI)
  message(
    FATAL_ERROR
      "MPI is required for hypre support. Set SUNDIALS_ENABLE_MPI to ON.")
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

find_package(HYPRE REQUIRED)

message(STATUS "HYPRE_LIBRARIES:   ${HYPRE_LIBRARIES}")
message(STATUS "HYPRE_INCLUDE_DIR: ${HYPRE_INCLUDE_DIR}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

if(SUNDIALS_ENABLE_HYPRE_CHECKS)

  message(CHECK_START "Testing hypre")

  # Create the test directory
  set(TEST_DIR ${PROJECT_BINARY_DIR}/HYPRE_TEST)

  # Create a C source file
  file(
    WRITE ${TEST_DIR}/test.c
    "\#include \"HYPRE_parcsr_ls.h\"\n"
    "int main(void) {\n"
    "HYPRE_ParVector par_b;\n"
    "HYPRE_IJVector b;\n"
    "par_b = 0;\n"
    "b = 0;\n"
    "if (par_b != 0 || b != 0) return(1);\n"
    "else return(0);\n"
    "}\n")

  # Attempt to build and link the test executable, pass --debug-trycompile to
  # the cmake command to save build files for debugging
  try_compile(
    COMPILE_OK ${TEST_DIR}
    ${TEST_DIR}/test.c
    LINK_LIBRARIES SUNDIALS::HYPRE MPI::MPI_C
    OUTPUT_VARIABLE COMPILE_OUTPUT)

  # Process test result
  if(COMPILE_OK)
    message(CHECK_PASS "success")
  else()
    message(CHECK_FAIL "failed")
    file(WRITE ${TEST_DIR}/compile.out "${COMPILE_OUTPUT}")
    message(
      FATAL_ERROR
        "Could not compile hypre test. Check output in ${TEST_DIR}/compile.out")
  endif()

else()
  message(STATUS "Skipped hypre tests.")
endif()
