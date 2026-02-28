# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ UMBC
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
# ARKode C++ parallel unit_tests
# ---------------------------------------------------------------

# List of test tuples of the form "name\;tasks\;args"
if(NOT SUNDIALS_PRECISION MATCHES "SINGLE")
  set(unit_tests "ark_test_heat2D_mri.cpp\;2\;0"
                 "ark_test_heat2D_mri.cpp\;4\;1")
endif()

# Add the build and install targets for each test
foreach(test_tuple ${unit_tests})

  # parse the test tuple
  list(GET test_tuple 0 test)
  list(GET test_tuple 1 number_of_tasks)
  list(GET test_tuple 2 test_args)

  # Extract the file name without extension
  get_filename_component(test_target ${test} NAME_WE)

  # check if this test has already been added, only need to add test source
  # files once for testing with different inputs
  if(NOT TARGET ${test_target})

    # test source files
    sundials_add_executable(${test_target} ${test})

    set_target_properties(${test_target} PROPERTIES FOLDER "unit_tests")

    # include location of public and private header files
    target_include_directories(
      ${test_target}
      PRIVATE $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
              ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src)

    # libraries to link against
    target_link_libraries(${test_target} MPI::MPI_CXX sundials_arkode
                          sundials_nvecparallel ${EXE_EXTRA_LINK_LIBS})

  endif()

  # check if test args are provided and set the test name
  if("${test_args}" STREQUAL "")
    set(test_name ${test_target})
  else()
    string(REGEX REPLACE " " "_" test_name "${test_target}_${test_args}")
  endif()

  # add test to regression tests
  sundials_add_test(
    ${test_name} ${test_target}
    TEST_ARGS ${test_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    MPI_NPROCS ${number_of_tasks})

endforeach()

message(STATUS "Added ARKODE CXX parallel unit tests")
