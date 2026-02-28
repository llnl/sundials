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
# Module to find and setup CALIPER correctly.
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

if(NOT DEFINED SUNDIALS_CALIPER_INCLUDED)
  set(SUNDIALS_CALIPER_INCLUDED)
else()
  return()
endif()

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

if(CALIPER_FOUND AND (NOT CALIPER_WORKS))
  # Do any checks which don't require compilation first.

  # Create the CALIPER_TEST directory
  set(CALIPER_TEST_DIR ${PROJECT_BINARY_DIR}/CALIPER_TEST)
  file(MAKE_DIRECTORY ${CALIPER_TEST_DIR})

  # If C++ is enabled, we build the example as a C++ code as a workaround for
  # what appears to be a bug in try_compile. If we dont do this, then
  # try_compile throws an error "No known features for CXX compiler".
  if(CXX_FOUND)
    set(_ext cpp)
  else()
    set(_ext c)
  endif()

  # Create a C source file
  file(
    WRITE ${CALIPER_TEST_DIR}/ltest.${_ext}
    "\#include <caliper/cali.h>\n"
    "int main(void)\n"
    "{\n"
    "  CALI_MARK_FUNCTION_BEGIN;\n"
    "  CALI_MARK_FUNCTION_END;\n"
    "  return 0;\n"
    "}\n")

  # Attempt to build and link executable with caliper
  try_compile(
    COMPILE_OK ${CALIPER_TEST_DIR}
    ${CALIPER_TEST_DIR}/ltest.${_ext}
    LINK_LIBRARIES caliper
    OUTPUT_VARIABLE COMPILE_OUTPUT)

  # Process test result
  if(COMPILE_OK)
    message(STATUS "Checking if CALIPER works with SUNDIALS... OK")
    set(CALIPER_WORKS
        TRUE
        CACHE BOOL "CALIPER works with SUNDIALS as configured" FORCE)
  else()
    message(STATUS "Checking if CALIPER works with SUNDIALS... FAILED")
    message(STATUS "Check output: ")
    message("${COMPILE_OUTPUT}")
    message(FATAL_ERROR "SUNDIALS interface to CALIPER is not functional.")
  endif()

elseif(CALIPER_FOUND AND CALIPER_WORKS)
  message(STATUS "Skipped CALIPER tests, assuming CALIPER works with SUNDIALS.")
endif()
