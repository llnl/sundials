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
set(_arkode_libs sundials_arkode sundials_nvecserial
                 sundials_sundomeigestpower ${EXE_EXTRA_LINK_LIBS})

sundials_add_example(
  ark_analytic ark_analytic.c
  INSTALL_FILES ark_analytic_arkode.scalar_tolerances_1e-6_1e-8_arkode.table_names_ARKODE_ESDIRK547L2SA_7_4_5_ARKODE_ERK_NONE.out
  LINK_LIBRARIES PRIVATE ${_arkode_libs})
sundials_add_test(
  ark_analytic_arkode.scalar_tolerances_1e-6_1e-8_arkode.table_names_ARKODE_ESDIRK547L2SA_7_4_5_ARKODE_ERK_NONE
  ark_analytic
  TEST_ARGS arkode.scalar_tolerances 1e-6 1e-8 arkode.table_names
             ARKODE_ESDIRK547L2SA_7_4_5 ARKODE_ERK_NONE
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE ark_analytic_arkode.scalar_tolerances_1e-6_1e-8_arkode.table_names_ARKODE_ESDIRK547L2SA_7_4_5_ARKODE_ERK_NONE.out
  EXAMPLE_TYPE develop)

foreach(_example ark_analytic_lsrk ark_analytic_lsrk_varjac
                 ark_analytic_lsrk_domeigest ark_analytic_mels
                 ark_analytic_nonlin ark_analytic_ssprk)
  if(_example STREQUAL "ark_analytic_lsrk_domeigest")
    sundials_add_example(
      ${_example} ${_example}.c EXAMPLE_TYPE develop
      INSTALL_FILES ark_analytic_lsrk_domeigest_arkid.dom_eig_est_init_preprocess_iters_1_sundomeigestimator.max_iters_1.out
      LINK_LIBRARIES PRIVATE ${_arkode_libs})
  else()
    sundials_add_example(${_example} ${_example}.c EXAMPLE_TYPE develop
                          LINK_LIBRARIES PRIVATE ${_arkode_libs})
  endif()
endforeach()

sundials_add_test(
  ark_analytic_lsrk_domeigest_arkid.dom_eig_est_init_preprocess_iters_1_sundomeigestimator.max_iters_1
  ark_analytic_lsrk_domeigest
  TEST_ARGS arkid.dom_eig_est_init_preprocess_iters 1 sundomeigestimator.max_iters 1
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE ark_analytic_lsrk_domeigest_arkid.dom_eig_est_init_preprocess_iters_1_sundomeigestimator.max_iters_1.out
  EXAMPLE_TYPE develop)

sundials_add_example(ark_analytic_partitioned ark_analytic_partitioned.c
                      TEST_NAME ark_analytic_partitioned_forcing
                      TEST_ARGS forcing EXAMPLE_TYPE develop
                      INSTALL_FILES ark_analytic_partitioned_splitting.out
                                     ark_analytic_partitioned_splitting_ARKODE_SPLITTING_BEST_2_2_2.out
                                     ark_analytic_partitioned_splitting_ARKODE_SPLITTING_RUTH_3_3_2.out
                                     ark_analytic_partitioned_splitting_ARKODE_SPLITTING_YOSHIDA_8_6_2.out
                      LINK_LIBRARIES PRIVATE ${_arkode_libs})
foreach(_case splitting
              "splitting ARKODE_SPLITTING_BEST_2_2_2"
              "splitting ARKODE_SPLITTING_RUTH_3_3_2"
              "splitting ARKODE_SPLITTING_YOSHIDA_8_6_2")
  string(REPLACE " " "_" _case_name "${_case}")
  set(_case_type develop)
  if(_case MATCHES "RUTH|YOSHIDA")
    set(_case_type exclude-single)
  endif()
  sundials_add_test(ark_analytic_partitioned_${_case_name}
                    ark_analytic_partitioned TEST_ARGS ${_case}
                    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                    ANSWER_FILE ark_analytic_partitioned_${_case_name}.out
                    EXAMPLE_TYPE ${_case_type})
endforeach()
