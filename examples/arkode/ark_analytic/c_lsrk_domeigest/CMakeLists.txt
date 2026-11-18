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
  ark_analytic_lsrk_domeigest
  SOURCES ark_analytic_lsrk_domeigest.c
  EXAMPLE_TYPE develop
  INSTALL_FILES
    ark_analytic_lsrk_domeigest_arkid.dom_eig_est_init_preprocess_iters_1_sundomeigestimator.max_iters_1.out
  LINK_LIBRARIES PRIVATE sundials_arkode sundials_nvecserial
                 sundials_sundomeigestpower ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(
  ark_analytic_lsrk_domeigest_arkid.dom_eig_est_init_preprocess_iters_1_sundomeigestimator.max_iters_1
  ark_analytic_lsrk_domeigest
  TEST_ARGS arkid.dom_eig_est_init_preprocess_iters 1
            sundomeigestimator.max_iters 1
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE
    ark_analytic_lsrk_domeigest_arkid.dom_eig_est_init_preprocess_iters_1_sundomeigestimator.max_iters_1.out
  EXAMPLE_TYPE develop)
