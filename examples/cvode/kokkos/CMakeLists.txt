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

# Example lists are tuples "name\;args\;type" where the type is 'develop' for
# examples excluded from 'make test' in releases
set(examples_list "cv_bruss_batched_kokkos.cpp\;\;develop"
                  "cv_bruss_batched_kokkos_2D.cpp\;\;develop")

# Add the build targets for each example
foreach(example_tuple ${examples_list})
  foreach(backend ${KOKKOS_EXAMPLES_BACKENDS})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 example_args)
    list(GET example_tuple 2 example_type)

    # extract the file name without extension
    get_filename_component(example_target ${example} NAME_WE)
    set(example_target "${example_target}.${backend}")

    if(NOT TARGET ${example_target})

      # example source files
      sundials_add_executable(${example_target} ${example})

      # folder for IDEs
      set_target_properties(${example_target} PROPERTIES FOLDER "Examples")

      # which backend to use
      target_compile_definitions(${example_target} PRIVATE USE_${backend})

      # directories to include
      target_include_directories(
        ${example_target} PRIVATE "${PROJECT_SOURCE_DIR}/examples/utilities")

      # libraries to link against
      target_link_libraries(
        ${example_target}
        PRIVATE sundials_cvode sundials_nveckokkos
                sundials_sunmatrixkokkosdense sundials_sunlinsolkokkosdense
                ${EXE_EXTRA_LINK_LIBS})

      if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
        target_link_libraries(${example_target}
                              PRIVATE sundials_cvode_fused_stubs)
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
      EXTRA_ARGS --kokkos-disable-warnings
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${test_name}.out
      EXAMPLE_TYPE ${example_type})

  endforeach()
endforeach()

# Install the targets
if(SUNDIALS_EXAMPLES_ENABLE_INSTALL)

  if(SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS)
    set(_fused_stubs cvode_fused_stubs)
  endif()

  sundials_install_examples(
    cvode examples_list
    CMAKE_TEMPLATE cmakelists_CXX_ex.in
    SUNDIALS_COMPONENTS nveckokkos sunmatrixkokkosdense sunlinsolkokkosdense
    SUNDIALS_TARGETS cvode ${_fused_stubs}
    OTHER_TARGETS Kokkos::kokkos Kokkos::kokkoskernels
    EXTRA_FILES "${PROJECT_SOURCE_DIR}/examples/utilities/example_utilities.hpp"
    DESTINATION cvode/kokkos)
endif()
