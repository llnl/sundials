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
# Module to find and setup ADIAK correctly.
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

if(NOT DEFINED SUNDIALS_ADIAK_INCLUDED)
  set(SUNDIALS_ADIAK_INCLUDED)
else()
  return()
endif()

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

if(adiak_FOUND AND (NOT adiak_WORKS))
  # Do any checks which don't require compilation first.

  # Create the adiak_TEST directory
  set(adiak_TEST_DIR ${PROJECT_BINARY_DIR}/adiak_TEST)
  file(MAKE_DIRECTORY ${adiak_TEST_DIR})

  # Create a C source file
  file(
    WRITE ${adiak_TEST_DIR}/ltest.c
    "\#include <adiak.h>\n"
    "int main(void)\n"
    "{\n"
    "  adiak_init(NULL);\n"
    "  adiak_fini();\n"
    "  return 0;\n"
    "}\n")

  # To ensure we do not use stuff from the previous attempts, we must remove the
  # CMakeFiles directory.
  file(REMOVE_RECURSE ${adiak_TEST_DIR}/CMakeFiles)

  # Attempt to build and link the "ltest" executable
  try_compile(
    COMPILE_OK ${adiak_TEST_DIR}
    ${adiak_TEST_DIR}/ltest.c
    OUTPUT_VARIABLE COMPILE_OUTPUT
    LINK_LIBRARIES adiak::adiak ${CMAKE_DL_LIBS})

  # Process test result
  if(COMPILE_OK)
    message(STATUS "Checking if adiak works with SUNDIALS... OK")
    set(adiak_WORKS
        TRUE
        CACHE BOOL "adiak works with SUNDIALS as configured" FORCE)
  else()
    message(STATUS "Checking if adiak works with SUNDIALS... FAILED")
    message(STATUS "Check output: ")
    message("${COMPILE_OUTPUT}")
    message(FATAL_ERROR "SUNDIALS interface to adiak is not functional.")
  endif()

elseif(adiak_FOUND AND adiak_WORKS)
  message(STATUS "Skipped adiak tests, assuming adiak works with SUNDIALS.")
endif()
