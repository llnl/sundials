# -----------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
# Module to find and setup XBraid.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

include_guard(GLOBAL)

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# Using XBRAID requires building with MPI enabled
if(NOT SUNDIALS_ENABLE_MPI)
  message(
    FATAL_ERROR
      "MPI is required for XBraid support. Set SUNDIALS_ENABLE_MPI to ON.")
endif()

# XBraid does not support single or extended precision
if(SUNDIALS_PRECISION MATCHES "SINGLE" OR SUNDIALS_PRECISION MATCHES "EXTENDED")
  message(
    FATAL_ERROR "XBraid is not compatible with ${SUNDIALS_PRECISION} precision")
endif()

# XBraid does not support 64-bit index sizes
if(SUNDIALS_INDEX_SIZE MATCHES "64")
  message(
    FATAL_ERROR
      "XBraid is not compatible with ${SUNDIALS_INDEX_SIZE}-bit indices")
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

find_package(XBRAID REQUIRED)

message(STATUS "XBRAID_LIBRARIES: ${XBRAID_LIBS}")
message(STATUS "XBRAID_INCLUDES:  ${XBRAID_INCS}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

if(SUNDIALS_ENABLE_XBRAID_CHECKS)

  message(CHECK_START "Testing XBraid")

  # Create the XBRAID_TEST directory
  set(TEST_DIR ${PROJECT_BINARY_DIR}/XBRAID_TEST)

  # Create a C source file
  file(
    WRITE ${TEST_DIR}/test.c
    "\#include <stdlib.h>\n"
    "\#include \"braid.h\"\n"
    "int main(void) {\n"
    "braid_Int rand;\n"
    "rand = braid_Rand();\n"
    "if (rand < 0) return 1;\n"
    "return 0;\n"
    "}\n")

  # Attempt to build and link the test executable, pass --debug-trycompile to
  # the cmake command to save build files for debugging
  try_compile(
    COMPILE_OK ${TEST_DIR}
    ${TEST_DIR}/test.c
    LINK_LIBRARIES SUNDIALS::XBRAID MPI::MPI_C ${SUNDIALS_MATH_LIBRARY}
    OUTPUT_VARIABLE COMPILE_OUTPUT)

  # Process test result
  if(COMPILE_OK)
    message(CHECK_PASS "success")
  else()
    message(CHECK_FAIL "failed")
    file(WRITE ${TEST_DIR}/compile.out "${COMPILE_OUTPUT}")
    message(
      FATAL_ERROR
        "Could not compile XBraid test. Check output in ${TEST_DIR}/compile.out"
    )
  endif()

else()
  message(STATUS "Skipped XBraid checks.")
endif()
