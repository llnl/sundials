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
# -----------------------------------------------------------------
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  cmake_minimum_required(VERSION 3.18)
  project(cvode_example LANGUAGES NONE)
endif()
include(${CMAKE_CURRENT_LIST_DIR}/../../SundialsCVODEExample.cmake)
# CVODE Ginkgo batched Brusselator example
foreach(backend REF OMP CUDA HIP SYCL)
  if(NOT SUNDIALS_GINKGO_BACKENDS MATCHES "${backend}")
    continue()
  endif()

  set(float_precision default)
  if(backend MATCHES "CUDA")
    set(vector nveccuda)
    set(language CUDA)
    set(float_precision 4)
  elseif(backend MATCHES "HIP")
    set(vector nvechip)
    set(language CXX)
  elseif(backend MATCHES "SYCL")
    set(vector nvecsycl)
    set(language CXX)
  elseif(backend MATCHES "OMP")
    set(vector nvecopenmp)
    set(language CXX)
  else()
    set(vector nvecserial)
    set(language CXX)
  endif()

  set_source_files_properties(cv_bruss_batched_ginkgo.cpp
                              PROPERTIES LANGUAGE ${language})
  set(example_target cv_bruss_batched_ginkgo.${backend})
  set(include_args)
  if(EXISTS "${PROJECT_SOURCE_DIR}/examples/utilities")
    set(include_args INCLUDE_DIRECTORIES PRIVATE
                     ${PROJECT_SOURCE_DIR}/examples/utilities)
  endif()
  sundials_add_executable(
    ${example_target}
    cv_bruss_batched_ginkgo.cpp
    LINK_LIBRARIES
    PRIVATE
    ${CVODE_C_LIBS}
    sundials_${vector}
    Ginkgo::ginkgo
    ${EXE_EXTRA_LINK_LIBS}
    ${include_args}
    COMPILE_DEFINITIONS
    PRIVATE
    USE_${backend}
    PROPERTIES
    FOLDER
    "Examples")
  sundials_add_test(
    ${example_target} ${example_target}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example_target}.out
    EXAMPLE_TYPE develop
    FLOAT_PRECISION ${float_precision})
endforeach()

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  file(GLOB example_out cv_bruss_batched_ginkgo*.out)
  install(
    FILES cv_bruss_batched_ginkgo.cpp cv_bruss_batched_ginkgo.hpp ${example_out}
          "${PROJECT_SOURCE_DIR}/examples/utilities/example_utilities.hpp"
    DESTINATION
      ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cv_brusselator/cpp-ginkgo)
endif()
