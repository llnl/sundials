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
#   sundials_install_examples_ginkgo(EXAMPLES_VAR
#     [TARGETS targets]
#     [BACKENDS backends]
#     [UNIT_TEST]
#   )
#
# adds a build target for each example tuple in EXAMPLES_VAR.
#
# The TARGETS option is a list of CMake targets provided to
# target_link_libraries.
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
  set(multiValueArgs TARGETS BACKENDS)

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

      if(NOT TARGET ${example_target})

        # create target
        add_executable(${example_target} ${example})

        # folder for IDEs
        set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

        # which backend to use
        target_compile_definitions(${example_target} PRIVATE USE_${backend})

        # directories to include
        target_include_directories(
          ${example_target} PRIVATE "${PROJECT_SOURCE_DIR}/examples/utilities")

        # libraries to link against
        target_link_libraries(
          ${example_target} PRIVATE ${arg_TARGETS} sundials_${vector}
                                    Ginkgo::ginkgo ${EXTRA_LINK_LIBS})

      endif()

      # check if example args are provided and set the test name
      if("${example_args}" STREQUAL "")
        set(test_name ${example_target})
      else()
        string(REGEX REPLACE " " "_" test_name
                             ${example_target}_${example_args})
      endif()

      # add example to regression tests
      if(${arg_UNIT_TEST})
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
