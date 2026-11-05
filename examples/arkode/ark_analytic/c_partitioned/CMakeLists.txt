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
  ark_analytic_partitioned
  SOURCES ark_analytic_partitioned.c
  TEST_NAME ark_analytic_partitioned_forcing
  TEST_ARGS forcing
  EXAMPLE_TYPE develop
  INSTALL_FILES
    ark_analytic_partitioned_splitting.out
    ark_analytic_partitioned_splitting_ARKODE_SPLITTING_BEST_2_2_2.out
    ark_analytic_partitioned_splitting_ARKODE_SPLITTING_RUTH_3_3_2.out
    ark_analytic_partitioned_splitting_ARKODE_SPLITTING_YOSHIDA_8_6_2.out
  LINK_LIBRARIES PRIVATE sundials_arkode sundials_nvecserial
                 sundials_sundomeigestpower ${EXE_EXTRA_LINK_LIBS})
foreach(
  _case
  splitting "splitting ARKODE_SPLITTING_BEST_2_2_2"
  "splitting ARKODE_SPLITTING_RUTH_3_3_2"
  "splitting ARKODE_SPLITTING_YOSHIDA_8_6_2")
  string(REPLACE " " "_" _case_name "${_case}")
  set(_case_type develop)
  if(_case MATCHES "RUTH|YOSHIDA")
    set(_case_type exclude-single)
  endif()
  sundials_add_test(
    ark_analytic_partitioned_${_case_name} ark_analytic_partitioned
    TEST_ARGS ${_case}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ark_analytic_partitioned_${_case_name}.out
    EXAMPLE_TYPE ${_case_type})
endforeach()
