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
sundials_add_example(
  cv_analytic_sys_dns_f2003
  cv_analytic_sys_dns_f2003.f90
  ADD_TEST
  ANSWER_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE
  cv_analytic_sys_dns_f2003.out
  EXAMPLE_TYPE
  develop
  INSTALL_FILES
  cv_analytic_sys_dns_f2003.out
  LINK_LIBRARIES
  ${CVODE_C_LIBS}
  sundials_fcvode_mod
  ${EXE_EXTRA_LINK_LIBS}
  PROPERTIES
  FOLDER
  "Examples"
  Fortran_MODULE_DIRECTORY
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cv_analytic_sys_dns_f2003.dir)

sundials_add_example(
  cv_analytic_sys_dns_jac_f2003
  cv_analytic_sys_dns_jac_f2003.f90
  ADD_TEST
  ANSWER_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE
  cv_analytic_sys_dns_jac_f2003.out
  EXAMPLE_TYPE
  develop
  INSTALL_FILES
  cv_analytic_sys_dns_jac_f2003.out
  LINK_LIBRARIES
  ${CVODE_C_LIBS}
  sundials_fcvode_mod
  ${EXE_EXTRA_LINK_LIBS}
  PROPERTIES
  FOLDER
  "Examples"
  Fortran_MODULE_DIRECTORY
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cv_analytic_sys_dns_jac_f2003.dir)
