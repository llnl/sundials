# -----------------------------------------------------------------------------
# Programmer(s): Yu Pan @ LLNL
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
# Module to find and setup ADIAK.
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

find_package(adiak REQUIRED)
message(STATUS "ADIAK_LIBRARIES:   ${adiak_LIBRARIES}")
message(STATUS "ADIAK_INCLUDE_DIR: ${adiak_INCLUDE_DIR}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

if(SUNDIALS_ENABLE_ADIAK_CHECKS)

  message(CHECK_START "Testing Adiak")

  # Create the adiak_TEST directory
  set(TEST_DIR ${PROJECT_BINARY_DIR}/ADIAK_TEST)

  # Create a C source file
  file(
    WRITE ${TEST_DIR}/test.c
    "\#include <adiak.h>\n"
    "int main(void)\n"
    "{\n"
    "  adiak_init(NULL);\n"
    "  adiak_fini();\n"
    "  return 0;\n"
    "}\n")

  # Attempt to build and link the test executable, pass --debug-trycompile to
  # the cmake command to save build files for debugging
  try_compile(
    COMPILE_OK ${TEST_DIR}
    ${TEST_DIR}/test.c
    LINK_LIBRARIES adiak::adiak ${CMAKE_DL_LIBS}
    OUTPUT_VARIABLE COMPILE_OUTPUT)

  # Check the result
  if(COMPILE_OK)
    message(CHECK_PASS "success")
  else()
    message(CHECK_FAIL "failed")
    file(WRITE ${TEST_DIR}/compile.out "${COMPILE_OUTPUT}")
    message(
      FATAL_ERROR
        "Could not compile Adiak test. Check output in ${TEST_DIR}/compile.out")
  endif()

else()
  message(STATUS "Skipped Adiak checks.")
endif()
