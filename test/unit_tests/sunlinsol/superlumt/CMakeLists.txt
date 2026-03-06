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
# CMakeLists.txt file for sunlinsol SuperLU_MT examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is 'develop' for
# examples excluded from 'make test' in releases

# Examples using the SuperLU_MT linear solver Note: threaded tests are excluded
# because of a potential bug in SuperLU_MT that causes random test failures
set(sunlinsol_superlumt_examples
    "test_sunlinsol_superlumt\;300 0 1 0\;"
    "test_sunlinsol_superlumt\;300 1 1 0\;"
    "test_sunlinsol_superlumt\;1000 0 3 0\;exclude"
    "test_sunlinsol_superlumt\;1000 1 3 0\;exclude")

# Dependencies for sunlinsol examples
set(sunlinsol_superlumt_dependencies test_sunlinsol)

# Add source directory to include directories
include_directories(. ..)

# Add the build and install targets for each example
foreach(example_tuple ${sunlinsol_superlumt_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # check if this example has already been added, only need to add example
  # source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    sundials_add_executable(${example} ${example}.c ../test_sunlinsol.c)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(
      ${example} sundials_nvecserial sundials_sunmatrixdense
      sundials_sunmatrixsparse sundials_sunlinsolsuperlumt
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

endforeach(example_tuple ${sunlinsol_superlumt_examples})
