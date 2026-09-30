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
sundials_add_executable(
  cv_roberts_dns_f2003
  cv_roberts_dns_f2003.f90
  LINK_LIBRARIES
  ${CVODE_C_LIBS}
  sundials_fcvode_mod
  ${EXE_EXTRA_LINK_LIBS}
  PROPERTIES
  FOLDER
  "Examples"
  Fortran_MODULE_DIRECTORY
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cv_roberts_dns_f2003.dir)
sundials_add_test(
  cv_roberts_dns_f2003 cv_roberts_dns_f2003
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE cv_roberts_dns_f2003.out
  EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(
    FILES cv_roberts_dns_f2003.f90 cv_roberts_dns_f2003.out
    DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cv_roberts/fortran)
endif()

sundials_add_executable(
  cv_roberts_dns_constraints_f2003
  cv_roberts_dns_constraints_f2003.f90
  LINK_LIBRARIES
  ${CVODE_C_LIBS}
  sundials_fcvode_mod
  ${EXE_EXTRA_LINK_LIBS}
  PROPERTIES
  FOLDER
  "Examples"
  Fortran_MODULE_DIRECTORY
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cv_roberts_dns_constraints_f2003.dir)
sundials_add_test(
  cv_roberts_dns_constraints_f2003 cv_roberts_dns_constraints_f2003
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE cv_roberts_dns_constraints_f2003.out
  EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(
    FILES cv_roberts_dns_constraints_f2003.f90
          cv_roberts_dns_constraints_f2003.out
    DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cv_roberts/fortran)
endif()
