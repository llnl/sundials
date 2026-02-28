# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
# KINSOL C++ serial unit_tests
# ---------------------------------------------------------------

# List of test tuples of the form "name\;args"
set(unit_tests "kin_test_getjac.cpp\;")

# Add the build and install targets for each test
foreach(test_tuple ${unit_tests})

  # Parse the test tuple
  list(GET test_tuple 0 test)
  list(GET test_tuple 1 test_args)

  # Extract the file name without extension
  get_filename_component(test_target ${test} NAME_WE)

  # Check if this test has already been added, only need to add test source
  # files once for testing with different inputs
  if(NOT TARGET ${test_target})

    # Test source files
    sundials_add_executable(${test_target} ${test})

    # Folder to organize targets in an IDE
    set_target_properties(${test_target} PROPERTIES FOLDER "unit_tests")

    # Include location of public and private header files
    target_include_directories(
      ${test_target}
      PRIVATE $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
              ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src)

    # Libraries to link against
    target_link_libraries(${test_target} sundials_kinsol sundials_nvecserial
                          ${EXE_EXTRA_LINK_LIBS})

  endif()

  # Check if test args are provided and set the test name
  if("${test_args}" STREQUAL "")
    set(test_name ${test_target})
  else()
    string(REPLACE " " "_" test_name "${test_target}_${test_args}")
  endif()

  # add test to regression tests
  sundials_add_test(
    ${test_name} ${test_target}
    TEST_ARGS ${test_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE "develop")

endforeach()

message(STATUS "Added KINSOL CXX serial units tests")
