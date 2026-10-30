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
set(_kin_analytic_outputs
    kin_analytic_fp_--damping_fn.out
    kin_analytic_fp_--damping_fp_0.5.out
    kin_analytic_fp_--m_aa_2.out
    kin_analytic_fp_--m_aa_2_--damping_aa_0.5.out
    kin_analytic_fp_--m_aa_2_--damping_fn.out
    kin_analytic_fp_--m_aa_2_--delay_aa_2.out
    kin_analytic_fp_--m_aa_2_--orth_aa_1.out
    kin_analytic_fp_--m_aa_2_--orth_aa_2.out
    kin_analytic_fp_--m_aa_2_--orth_aa_3.out
    kin_analytic_fp_--m_aa_3_--depth_fn.out)
sundials_add_example(
  kin_analytic_fp SOURCES kin_analytic_fp.c NO_TEST
  INSTALL_FILES kin_analytic_fp.out ${_kin_analytic_outputs}
  LINK_LIBRARIES PRIVATE ${KINSOL_C_LIBS} sundials_nvecserial
                 ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(
  kin_analytic_fp kin_analytic_fp
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE kin_analytic_fp.out)

set(_kin_analytic_tests
    "--damping_fp 0.5"
    "--damping_fn"
    "--m_aa 2"
    "--m_aa 2 --delay_aa 2"
    "--m_aa 2 --damping_aa 0.5"
    "--m_aa 2 --damping_fn"
    "--m_aa 3 --depth_fn"
    "--m_aa 2 --orth_aa 1"
    "--m_aa 2 --orth_aa 2"
    "--m_aa 2 --orth_aa 3")
foreach(test_args ${_kin_analytic_tests})
  string(REPLACE " " "_" test_suffix "${test_args}")
  set(test_options)
  if(NOT test_args STREQUAL "--m_aa 2"
     AND NOT test_args STREQUAL "--m_aa 2 --delay_aa 2"
     AND NOT test_args MATCHES "--orth_aa")
    list(APPEND test_options EXAMPLE_TYPE develop)
  endif()
  sundials_add_test(
    kin_analytic_fp_${test_suffix} kin_analytic_fp
    TEST_ARGS ${test_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE kin_analytic_fp_${test_suffix}.out ${test_options})
endforeach()
