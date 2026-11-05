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
  ark_analytic
  SOURCES ark_analytic.c
  NO_TEST
  INSTALL_FILES
    ark_analytic.out
    ark_analytic_arkode.scalar_tolerances_1e-6_1e-8_arkode.table_names_ARKODE_ESDIRK547L2SA_7_4_5_ARKODE_ERK_NONE.out
  LINK_LIBRARIES PRIVATE sundials_arkode sundials_nvecserial
                 sundials_sundomeigestpower ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(
  ark_analytic ark_analytic
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE ark_analytic.out)
sundials_add_test(
  ark_analytic_arkode.scalar_tolerances_1e-6_1e-8_arkode.table_names_ARKODE_ESDIRK547L2SA_7_4_5_ARKODE_ERK_NONE
  ark_analytic
  TEST_ARGS arkode.scalar_tolerances 1e-6 1e-8 arkode.table_names
            ARKODE_ESDIRK547L2SA_7_4_5 ARKODE_ERK_NONE
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE
    ark_analytic_arkode.scalar_tolerances_1e-6_1e-8_arkode.table_names_ARKODE_ESDIRK547L2SA_7_4_5_ARKODE_ERK_NONE.out
)
