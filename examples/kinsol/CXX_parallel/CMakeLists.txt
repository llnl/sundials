# ------------------------------------------------------------------------------
# Programmer(s): Shelby Lockhart @ UIUC/LLNL
# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------
# CMakeLists.txt file for Kinsol MPI examples
# ------------------------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;tasks\;type" where the type is
# develop for examples excluded from 'make test' in releases.

# Examples to build
set(KINSOL_examples "kin_heat2D_nonlin_p.cpp\;--np 2 2\;1\;4\;exclude-single"
                    "kin_em_p.cpp\;\;1\;2\;exclude-single")

# Add the build target for each example
foreach(example_tuple ${KINSOL_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 number_of_nodes)
  list(GET example_tuple 3 number_of_tasks)
  list(GET example_tuple 4 example_type)

  # extract the file name without extension
  get_filename_component(example_target ${example} NAME_WE)

  if(NOT TARGET ${example_target})

    # add example source files
    sundials_add_executable(${example_target} ${example})

    set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(
      ${example_target}
      PRIVATE sundials_kinsol sundials_nvecparallel sundials_nvecmpiplusx
              sundials_nvecserial MPI::MPI_CXX ${EXE_EXTRA_LINK_LIBS})

  endif()

  # set the test name
  if(example_args)
    string(REGEX REPLACE " " "_" test_name ${example_target}_${example_args})
  else()
    set(test_name ${example_target})
  endif()

  # add regression test
  sundials_add_test(
    ${test_name} ${example_target}
    TEST_ARGS ${example_args}
    MPI_NPROCS ${number_of_tasks}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE ${example_type})

endforeach()

# create Makefile and CMakeLists.txt for examples
if(SUNDIALS_EXAMPLES_ENABLE_INSTALL)

  sundials_install_examples(
    kinsol KINSOL_examples
    CMAKE_TEMPLATE cmakelists_CXX_MPI_ex.in
    MAKE_TEMPLATE makefile_parallel_CXX_ex.in
    SOLVER_LIBRARY sundials_kinsol
    SUNDIALS_TARGETS kinsol nvecparallel nvecmpiplusx nvecserial
    DESTINATION kinsol/CXX_parallel
    EXTRA_FILES README
    TEST_INSTALL CXX_parallel)

endif()
