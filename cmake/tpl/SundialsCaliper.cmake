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
# Module to find and setup CALIPER.
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

find_package(CALIPER PATHS "${CALIPER_DIR}" REQUIRED)

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

if(SUNDIALS_ENABLE_CALIPER_CHECKS)

  message(CHECK_START "Testing Caliper")

  # Create the test directory
  set(TEST_DIR ${PROJECT_BINARY_DIR}/CALIPER_TEST)

  # If C++ is enabled, we build the example as a C++ code as a workaround for
  # what appears to be a bug in try_compile. If we don't do this, then
  # try_compile throws an error "No known features for CXX compiler".
  if(CXX_FOUND)
    set(_ext cpp)
  else()
    set(_ext c)
  endif()

  # Create a source test file
  file(
    WRITE ${TEST_DIR}/test.${_ext}
    "\#include <caliper/cali.h>\n"
    "int main(void)\n"
    "{\n"
    "  CALI_MARK_FUNCTION_BEGIN;\n"
    "  CALI_MARK_FUNCTION_END;\n"
    "  return 0;\n"
    "}\n")

  # Attempt to build and link the test executable, pass --debug-trycompile to
  # the cmake command to save build files for debugging
  try_compile(
    COMPILE_OK ${TEST_DIR}
    ${TEST_DIR}/test.${_ext}
    LINK_LIBRARIES caliper
    OUTPUT_VARIABLE COMPILE_OUTPUT)

  # Check the result
  if(COMPILE_OK)
    message(CHECK_PASS "success")
  else()
    message(CHECK_FAIL "failed")
    file(WRITE ${TEST_DIR}/compile.out "${COMPILE_OUTPUT}")
    message(
      FATAL_ERROR
        "Could not compile Caliper test. Check output in ${TEST_DIR}/compile.out"
    )
  endif()

else()
  message(STATUS "Skipped Caliper checks.")
endif()
