# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------
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
# CMakeLists.txt file for slunrloc_sunmatrix examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;type" where the 'develop' for
# examples excluded from 'make test' in releases

# Examples using SUNDIALS SUNMatrix wrapper of the SuperLU SLU_NR_loc
# SuperMatrix.
set(sunmatrix_slunrloc_examples
    "test_sunmatrix_slunrloc.cpp\;100 1 1 0 0\;1\;\;"
    "test_sunmatrix_slunrloc.cpp\;400 2 1 0 0\;2\;\;"
    "test_sunmatrix_slunrloc.cpp\;600 3 1 0 0\;3\;\;"
    "test_sunmatrix_slunrloc.cpp\;1000 1 4 0 0\;4\;\;")

# Dependencies for sunmatrix examples
set(sunmatrix_slunrloc_dependencies test_sunmatrix)

# Add the build and install targets for each example
foreach(example_tuple ${sunmatrix_slunrloc_examples})
  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 number_of_nodes)
  list(GET example_tuple 3 example_type)

  # extract the file name without extension
  get_filename_component(example_target ${example} NAME_WE)

  # check if this example has already been added, only need to add example
  # source files once for testing with different inputs
  if(NOT TARGET ${example_target})
    # example source files
    sundials_add_executable(${example_target} ${example} ../test_sunmatrix.c)

    # folder to organize targets in an IDE
    set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

    target_include_directories(${example_target} PRIVATE . ..)

    # libraries to link against
    target_link_libraries(
      ${example_target}
      PRIVATE MPI::MPI_CXX sundials_sunmatrixslunrloc sundials_sunmatrixdense
              sundials_nvecparallel sundials_nvecserial)
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    if("${number_of_nodes}" STREQUAL "")
      set(test_name ${example_target})
    else()
      string(REGEX REPLACE " " "_" test_name
                           ${example_target}_${number_of_nodes})
    endif()
  else()
    if("${number_of_nodes}" STREQUAL "")
      string(REGEX REPLACE " " "_" test_name ${example_target}_${example_args})
    else()
      string(REGEX REPLACE " " "_" test_name
                           ${example_target}_${number_of_nodes}_${example_args})
    endif()
  endif()

  # add example to regression tests
  sundials_add_test(
    ${test_name} ${example_target}
    TEST_ARGS ${example_args}
    MPI_NPROCS ${number_of_nodes}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

endforeach(example_tuple ${sunmatrix_slunrloc_examples})
