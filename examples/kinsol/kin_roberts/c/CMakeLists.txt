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
  kin_roberts_fp SOURCES kin_roberts_fp.c
  INSTALL_FILES kin_roberts_fp_kinsol.m_aa_1.out
  LINK_LIBRARIES PRIVATE ${KINSOL_C_LIBS} sundials_nvecserial
                 ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(
  kin_roberts_fp_kinsol.m_aa_1 kin_roberts_fp
  TEST_ARGS kinsol.m_aa 1
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE kin_roberts_fp_kinsol.m_aa_1.out
  EXAMPLE_TYPE develop)
