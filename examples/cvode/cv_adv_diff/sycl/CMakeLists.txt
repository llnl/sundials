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
sundials_add_executable(cv_adv_diff_kry_sycl cv_adv_diff_kry_sycl.cpp)
set_target_properties(cv_adv_diff_kry_sycl PROPERTIES FOLDER "Examples")
target_link_libraries(
  cv_adv_diff_kry_sycl PRIVATE sundials_cvode sundials_nvecsycl
                               ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(
  cv_adv_diff_kry_sycl cv_adv_diff_kry_sycl
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE cv_adv_diff_kry_sycl.out
  EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cv_adv_diff_kry_sycl.cpp cv_adv_diff_kry_sycl.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cv_adv_diff/sycl)
endif()
