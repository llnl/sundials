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
  cvs_analytic_mels cvs_analytic_mels.c
  INSTALL_FILES cvs_analytic_mels_cvodes.max_order_3.out
  LINK_LIBRARIES PRIVATE ${CVODES_C_LIBS} sundials_nvecserial
                 ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(
  cvs_analytic_mels_cvodes.max_order_3 cvs_analytic_mels
  TEST_ARGS cvodes.max_order 3
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE cvs_analytic_mels_cvodes.max_order_3.out
  EXAMPLE_TYPE develop)
