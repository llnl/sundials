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
# CVODE serial C++ KPR example
sundials_add_executable(cv_kpr cv_kpr.cpp)
set_target_properties(cv_kpr PROPERTIES FOLDER "Examples")
if(EXISTS "${PROJECT_SOURCE_DIR}/examples/utilities")
  target_include_directories(cv_kpr
                             PRIVATE ${PROJECT_SOURCE_DIR}/examples/utilities)
else()
  target_include_directories(cv_kpr PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
endif()
target_link_libraries(cv_kpr PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                                     ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(
  cv_kpr cv_kpr
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE cv_kpr.out
  EXAMPLE_TYPE develop)
sundials_add_test(
  cv_kpr_cvode.init_step_0.01 cv_kpr
  TEST_ARGS cvode.init_step 0.01
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE cv_kpr_cvode.init_step_0.01.out
  EXAMPLE_TYPE develop)

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cv_kpr.cpp cv_kpr.hpp cv_kpr.out cv_kpr_cvode.init_step_0.01.out
                "${PROJECT_SOURCE_DIR}/examples/utilities/example_utilities.hpp"
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cv_kpr/cpp)
endif()
