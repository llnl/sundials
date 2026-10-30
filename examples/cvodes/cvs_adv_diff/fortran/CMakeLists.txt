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
  cvs_adv_diff_fsa_non_f2003 SOURCES cvs_adv_diff_fsa_non_f2003.f90
  TEST_NAME cvs_adv_diff_fsa_non_f2003_-sensi_stg_t
  TEST_ARGS -sensi stg t
  INSTALL_FILES cvs_adv_diff_fsa_non_f2003_-sensi_sim_t.out
  LINK_LIBRARIES PRIVATE ${CVODES_C_LIBS} sundials_fcvodes_mod
                 ${EXE_EXTRA_LINK_LIBS}
  PROPERTIES
    Fortran_MODULE_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/cvs_adv_diff_fsa_non_f2003.dir)
sundials_add_test(
  cvs_adv_diff_fsa_non_f2003_-sensi_sim_t cvs_adv_diff_fsa_non_f2003
  TEST_ARGS -sensi sim t
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE cvs_adv_diff_fsa_non_f2003_-sensi_sim_t.out
  EXAMPLE_TYPE develop)
