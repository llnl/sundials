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
# CMakeLists.txt file for ARKODE C++ manyvector examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is 'develop' for
# examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(ARKODE_examples "ark_sod_lsrk.cpp\;\;exclude-single")

# Header files to install
set(ARKODE_headers)

# Auxiliary files to install
set(ARKODE_extras plot_sod.py)

# Add the build and install targets for each example
foreach(example_tuple ${ARKODE_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # extract the file name without extension
  get_filename_component(example_target ${example} NAME_WE)

  # add example executable
  if(NOT TARGET ${example_target})
    sundials_add_executable(${example_target} ${example})

    set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

    # directories to include
    target_include_directories(${example_target}
                               PRIVATE ${PROJECT_SOURCE_DIR}/examples/utilities)

    # libraries to link against
    target_link_libraries(${example_target} sundials_arkode sundials_nvecserial
                          sundials_nvecmanyvector ${EXE_EXTRA_LINK_LIBS})

  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    set(test_name ${example_target})
  else()
    string(REGEX REPLACE " " "_" test_name ${example_target}_${example_args})
  endif()

  # add example to regression tests
  sundials_add_test(
    ${test_name} ${example_target}
    TEST_ARGS ${example_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE ${example_type})

endforeach(example_tuple ${ARKODE_examples})

if(SUNDIALS_EXAMPLES_ENABLE_INSTALL)

  sundials_install_examples(
    arkode ARKODE_examples
    CMAKE_TEMPLATE cmakelists_CXX_ex.in
    MAKE_TEMPLATE makefile_serial_CXX_ex.in
    SUNDIALS_TARGETS nvecserial nvecmanyvector arkode
    DESTINATION arkode/CXX_manyvector
    EXTRA_FILES ${ARKODE_extras} ${ARKODE_headers}
    TEST_INSTALL CXX_manyvector)

endif()
