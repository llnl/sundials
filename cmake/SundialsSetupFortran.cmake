# ---------------------------------------------------------------
# Programmer(s): Radu Serban, David Gardner, Cody J. Balos @ LLNL
# ---------------------------------------------------------------
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
# ---------------------------------------------------------------
# Module which enables Fortran and tests for support of necessary
# compiler features for the current SUNDIALS configuration.
# Will define the variables:
#   Fortran_FOUND - TRUE if a Fortran compiler is found
#   F2003_FOUND   - TRUE if the Fortran compiler supports the
#                   Fortran 2003 standard
# ---------------------------------------------------------------

# -----------------------------------------------------------------------------
# Enable Fortran
# -----------------------------------------------------------------------------
enable_language(Fortran)
set(Fortran_FOUND TRUE)

# Enable preprocessing Fortran code
set(CMAKE_Fortran_PREPROCESS ON)

# -----------------------------------------------------------------------------
# Check if Fortran 2003 is supported
# -----------------------------------------------------------------------------
if(BUILD_FORTRAN_MODULE_INTERFACE)
  if(NOT F2003_FOUND)
    message(STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports F2003")

    set(F2003Test_DIR ${PROJECT_BINARY_DIR}/F2003Test_DIR)
    file(MAKE_DIRECTORY ${F2003Test_DIR})

    # Create a CMakeLists.txt file
    file(
      WRITE ${F2003Test_DIR}/CMakeLists.txt
      "CMAKE_MINIMUM_REQUIRED(VERSION ${CMAKE_VERSION})\n"
      "PROJECT(ftest Fortran)\n"
      "SET(CMAKE_VERBOSE_MAKEFILE ON)\n"
      "SET(CMAKE_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\")\n"
      "SET(CMAKE_Fortran_COMPILER \"${CMAKE_Fortran_COMPILER}\")\n"
      "SET(CMAKE_Fortran_FLAGS \"${CMAKE_Fortran_FLAGS}\")\n"
      "SET(CMAKE_Fortran_FLAGS_RELEASE \"${CMAKE_Fortran_FLAGS_RELEASE}\")\n"
      "SET(CMAKE_Fortran_FLAGS_DEBUG \"${CMAKE_Fortran_FLAGS_DEBUG}\")\n"
      "SET(CMAKE_Fortran_FLAGS_RELWITHDEBUGINFO \"${CMAKE_Fortran_FLAGS_RELWITHDEBUGINFO}\")\n"
      "SET(CMAKE_Fortran_FLAGS_MINSIZE \"${CMAKE_Fortran_FLAGS_MINSIZE}\")\n"
      "ADD_EXECUTABLE(ftest ftest.f90)\n")

    # Create a Fortran source file which tries to use iso_c_binding
    file(WRITE ${F2003Test_DIR}/ftest.f90
         "program main\n" "use, intrinsic :: iso_c_binding\n"
         "end program main\n")

    # Attempt compile the executable
    try_compile(
      FTEST_OK ${F2003Test_DIR}
      ${F2003Test_DIR} ftest
      OUTPUT_VARIABLE COMPILE_OUTPUT)

    # To ensure we do not use stuff from the previous attempts, we must remove
    # the CMakeFiles directory.
    file(REMOVE_RECURSE ${F2003Test_DIR}/CMakeFiles)

    if(FTEST_OK)
      message(
        STATUS
          "Checking whether ${CMAKE_Fortran_COMPILER} supports F2003 -- yes")
      set(F2003_FOUND
          TRUE
          CACHE BOOL "${CMAKE_Fortran_COMPILER} supports F2003" FORCE)
    else()
      message(
        STATUS "Checking whether ${CMAKE_Fortran_COMPILER} supports F2003 -- no"
      )
      message(STATUS "Check output:")
      message("${COMPILE_OUTPUT}")
      message(
        FATAL_ERROR
          "BUILD_FORTRAN_MODULE_INTERFACE is set to ON, but the CMAKE_Fortran_COMPILER does not support F2003"
      )
    endif()
  else()
    message(
      STATUS
        "Skipped F2003 tests, assuming ${CMAKE_Fortran_COMPILER} supports the f2003 standard. To rerun the F2003 tests, set F2003_FOUND to FALSE."
    )
  endif()
endif()

# Put all F2003 modules into one build directory
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/fortran")
