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
# CMakeLists.txt file for CVODE C++ serial examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;tasks\;type" where the type is
# develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(CVODE_examples
    "cv_heat2D.cpp\;\;\;\;exclude-single" "cv_kpr.cpp\;\;\;\;develop"
    "cv_kpr.cpp\;cvode.init_step 0.01\;\;\;develop")

# Add the build and install targets for each example
foreach(example_tuple ${CVODE_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 number_of_nodes)
  list(GET example_tuple 3 number_of_tasks)
  list(GET example_tuple 4 example_type)

  # extract the file name without extension
  get_filename_component(example_target ${example} NAME_WE)

  # check if this example has already been added, only need to add example
  # source files once for testing with different inputs
  if(NOT TARGET ${example_target})
    # example source files
    sundials_add_executable(${example_target} ${example})

    # folder to organize targets in an IDE
    set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

    # directories to include
    target_include_directories(${example_target}
                               PRIVATE ${PROJECT_SOURCE_DIR}/examples/utilities)

    # libraries to link against
    target_link_libraries(${example_target} sundials_cvode sundials_nvecserial
                          ${EXE_EXTRA_LINK_LIBS})

    if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
      target_link_libraries(${example_target} sundials_cvode_fused_stubs)
    endif()
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

endforeach(example_tuple ${CVODE_examples})

# create Makefile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
    set(_fused_stubs cvode_fused_stubs)
  endif()

  sundials_install_examples(
    cvode CVODE_examples
    CMAKE_TEMPLATE cmakelists_CXX_ex.in
    MAKE_TEMPLATE makefile_serial_CXX_ex.in
    SUNDIALS_TARGETS cvode nvecserial ${_fused_stubs}
    DESTINATION cvode/CXX_serial
    EXTRA_FILES README plot_heat2D.py
                "${PROJECT_SOURCE_DIR}/examples/utilities/example_utilities.hpp"
    TEST_INSTALL CXX_serial)

endif(EXAMPLES_INSTALL)
