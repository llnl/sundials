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
# KINSOL C serial unit tests
# ---------------------------------------------------------------

# List of test tuples of the form "name\;args"
set(unit_tests "kin_test_getuserdata\;" "kin_test_reuse_fp\;0"
               "kin_test_reuse_fp\;1")

# Add the build and install targets for each test
foreach(test_tuple ${unit_tests})

  # parse the test tuple
  list(GET test_tuple 0 test)
  list(GET test_tuple 1 test_args)

  # check if this test has already been added, only need to add test source
  # files once for testing with different inputs
  if(NOT TARGET ${test})

    # test source files
    sundials_add_executable(${test} ${test}.c)

    set_target_properties(${test} PROPERTIES FOLDER "unit_tests")

    # include location of public and private header files
    target_include_directories(
      ${test} PRIVATE $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
                      ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src)

    # libraries to link against
    target_link_libraries(${test} sundials_kinsol sundials_nvecserial
                          ${EXE_EXTRA_LINK_LIBS})

  endif()

  # check if test args are provided and set the test name
  if("${test_args}" STREQUAL "")
    set(test_name ${test})
  else()
    string(REPLACE " " "_" test_name "${test}_${test_args}")
    string(REPLACE " " ";" test_args "${test_args}")
  endif()

  # add test to regression tests
  add_test(NAME ${test_name} COMMAND ${test} ${test_args})

endforeach()

message(STATUS "Added KINSOL C serial units tests")
