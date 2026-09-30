# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
# The macro:
#
#   sundials_add_examples_ginkgo(EXAMPLES_VAR
#     [TARGETS targets]
#     [BACKENDS backends]
#     [INSTALL_FILES files]
#     [UNIT_TEST]
#   )
#
# adds a build target for each example tuple in EXAMPLES_VAR.
#
# The TARGETS option is a list of CMake targets provided to
# target_link_libraries. INSTALL_FILES are passed to sundials_add_example for
# installation alongside each example source and answer file.
#
# The BACKENDS is a list of Ginkgo backends compatible with the examples in
# EXAMPLES_VAR.
#
# When the UNIT_TEST option is provided sundials_add_test is called with NODIFF
# so the example return value determines pass/fail. Otherwise, the ANSWER_DIR
# and ANSWER_FILE are used to set an output file for comparison.
# ------------------------------------------------------------------------------

# Add the build targets for each CVODE example
macro(sundials_add_examples_ginkgo EXAMPLES_VAR)

  set(options UNIT_TEST)
  set(oneValueArgs)
  set(multiValueArgs TARGETS BACKENDS INSTALL_FILES)

  # Parse keyword arguments and options
  cmake_parse_arguments(arg "${options}" "${oneValueArgs}" "${multiValueArgs}"
                        ${ARGN})

  foreach(example_tuple ${${EXAMPLES_VAR}})
    foreach(backend ${arg_BACKENDS})

      # parse the example tuple
      list(GET example_tuple 0 example)
      list(GET example_tuple 1 example_args)
      list(GET example_tuple 2 example_type)

      if(NOT (SUNDIALS_GINKGO_BACKENDS MATCHES "${backend}"))
        continue()
      endif()

      set(float_precision "default")
      if(backend MATCHES "CUDA")
        set(vector nveccuda)
        set(float_precision "4")
      elseif(backend MATCHES "HIP")
        set(vector nvechip)
      elseif(backend MATCHES "SYCL")
        set(vector nvecsycl)
      elseif(backend MATCHES "OMP")
        set(vector nvecopenmp)
      elseif(backend MATCHES "REF")
        set(vector nvecserial)
      endif()

      if(backend MATCHES "CUDA")
        set_source_files_properties(${example} PROPERTIES LANGUAGE CUDA)
      else()
        set_source_files_properties(${example} PROPERTIES LANGUAGE CXX)
      endif()

      # extract the file name without extension
      get_filename_component(example_target ${example} NAME_WE)
      set(example_target "${example_target}.${backend}")

      # check if example args are provided and set the test name
      if("${example_args}" STREQUAL "")
        set(test_name ${example_target})
      else()
        string(REGEX REPLACE " " "_" test_name
                             ${example_target}_${example_args})
      endif()

      if(NOT TARGET ${example_target})
        set(test_args
            ADD_TEST
            TEST_NAME
            ${test_name}
            EXAMPLE_TYPE
            ${example_type}
            TEST_ARGS
            ${example_args})
        set(install_files ${arg_INSTALL_FILES})
        if(${arg_UNIT_TEST})
          list(APPEND test_args NODIFF)
        else()
          list(
            APPEND
            test_args
            ANSWER_DIR
            ${CMAKE_CURRENT_SOURCE_DIR}
            ANSWER_FILE
            ${test_name}.out
            FLOAT_PRECISION
            ${float_precision})
          list(APPEND install_files ${test_name}.out)
        endif()
        if(EXISTS "${PROJECT_SOURCE_DIR}/examples/utilities")
          set(example_utilities_dir "${PROJECT_SOURCE_DIR}/examples/utilities")
        else()
          set(example_utilities_dir "${CMAKE_CURRENT_SOURCE_DIR}")
        endif()
        sundials_add_example(
          ${example_target} ${example} ${test_args} INSTALL
          INSTALL_FILES ${install_files}
          LINK_LIBRARIES PRIVATE ${arg_TARGETS} sundials_${vector}
                         Ginkgo::ginkgo ${EXTRA_LINK_LIBS}
          INCLUDE_DIRECTORIES PRIVATE "${example_utilities_dir}"
          COMPILE_DEFINITIONS PRIVATE USE_${backend}
          PROPERTIES FOLDER "Examples")
      elseif(${arg_UNIT_TEST})
        sundials_add_test(
          ${test_name} ${example_target}
          EXAMPLE_TYPE ${example_type}
          TEST_ARGS ${example_args}
          NODIFF)
      else()
        sundials_add_test(
          ${test_name} ${example_target}
          EXAMPLE_TYPE ${example_type}
          TEST_ARGS ${example_args}
          ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
          ANSWER_FILE ${test_name}.out
          FLOAT_PRECISION ${float_precision})
      endif()

    endforeach()
  endforeach()

endmacro()
