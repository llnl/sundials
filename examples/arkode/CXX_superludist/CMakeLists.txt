# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# -----------------------------------------------------------------
# CMakeLists.txt file for ARKODE C++ SuperLU_DIST examples
# -----------------------------------------------------------------

# Example lists are tuples "name\;nodes\;tasks\;type" where the type is develop
# for examples excluded from 'make test' in releases

list(APPEND ARKODE_examples "ark_brusselator1D_FEM_sludist.cpp\;1\;1\;develop")

# Auxiliary files to install
set(ARKODE_extras)

# Add the build and install targets for each example
foreach(example_tuple ${ARKODE_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 number_of_nodes)
  list(GET example_tuple 2 number_of_tasks)
  list(GET example_tuple 3 example_type)

  # extract the file name without extension
  get_filename_component(example_target ${example} NAME_WE)

  # example source files
  sundials_add_executable(${example_target} ${example})

  set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

  # add example to regression tests
  sundials_add_test(
    ${example_target} ${example_target}
    MPI_NPROCS ${number_of_tasks}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example_target}.out
    EXAMPLE_TYPE ${example_type})

  # libraries to link against
  target_link_libraries(${example_target} PRIVATE MPI::MPI_CXX sundials_arkode
                                                  sundials_sunlinsolsuperludist)

endforeach(example_tuple ${ARKODE_examples})

if(SUPERLUDIST_CUDA)
  set(_ex_lang CUDA)
elseif(SUPERLUDIST_ROCM)
  set(_ex_lang HIP)
else()
  set(_ex_lang CXX)
endif()

# create Makefile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  sundials_install_examples(
    arkode_superludist ARKODE_examples
    EXTRA_FILES README ${ARKODE_extras}
    CMAKE_TEMPLATE cmakelists_${_ex_lang}_MPI_ex.in
    SUNDIALS_TARGETS arkode sunlinsolsuperludist
    DESTINATION arkode/CXX_superludist)

  # add test_install target
  sundials_add_test_install(arkode CXX_superludist)

endif()
