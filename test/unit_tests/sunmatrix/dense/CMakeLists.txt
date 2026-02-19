# ---------------------------------------------------------------
# Programmer(s): Daniel Reynolds @ UMBC
#                David J. Gardner and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for dense sunmatrix examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is 'develop' for
# examples excluded from 'make test' in releases

# Examples using SUNDIALS dense matrix
set(sunmatrix_dense_examples
    "test_sunmatrix_dense\;100 100 0\;" "test_sunmatrix_dense\;200 1000 0\;"
    "test_sunmatrix_dense\;2000 100 0\;")

# Dependencies for sunmatrix examples
set(sunmatrix_dense_dependencies test_sunmatrix)

# If building F2003 tests
if(BUILD_FORTRAN_MODULE_INTERFACE)
  set(sunmatrix_dense_fortran_examples "test_fsunmatrix_dense_mod\;\;")
endif()

# Add source directory to include directories
include_directories(. ..)

# Add the build and install targets for each example
foreach(example_tuple ${sunmatrix_dense_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # check if this example has already been added, only need to add example
  # source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    sundials_add_executable(${example} ${example}.c ../test_sunmatrix.c)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example} sundials_nvecserial
                          sundials_sunmatrixdense ${EXE_EXTRA_LINK_LIBS})
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    set(test_name ${example})
  else()
    string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
  endif()

  # add example to regression tests
  sundials_add_test(
    ${test_name} ${example}
    TEST_ARGS ${example_args}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

endforeach(example_tuple ${sunmatrix_dense_examples})

# Add the build and install targets for each example
foreach(example_tuple ${sunmatrix_dense_fortran_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # check if this example has already been added, only need to add example
  # source files once for testing with different inputs
  if(NOT TARGET ${example})
    # build fortran modules into a unique directory to avoid naming collisions
    set(CMAKE_Fortran_MODULE_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${example}.dir)

    sundials_add_executable(
      ${example} ${example}.f90
      ${SUNDIALS_SOURCE_DIR}/test/unit_tests/utilities/test_utilities.f90
      ${SUNDIALS_SOURCE_DIR}/test/unit_tests/sunmatrix/test_sunmatrix.f90)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(
      ${example} sundials_nvecserial sundials_fnvecserial_mod
      sundials_sunmatrixdense sundials_fsunmatrixdense_mod
      ${EXE_EXTRA_LINK_LIBS})
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    set(test_name ${example})
  else()
    string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
  endif()

  # add example to regression tests
  sundials_add_test(
    ${test_name} ${example}
    TEST_ARGS ${example_args}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

endforeach(example_tuple ${sunmatrix_dense_fortran_examples})
