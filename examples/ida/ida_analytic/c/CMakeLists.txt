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
  idaAnalytic_mels idaAnalytic_mels.c
  EXAMPLE_TYPE exclude-single
  INSTALL_FILES idaAnalytic_mels_ida.scalar_tolerances_1e-3_1e-8.out
  LINK_LIBRARIES PRIVATE ${IDA_C_LIBS} sundials_nvecserial
                 ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(
  idaAnalytic_mels_ida.scalar_tolerances_1e-3_1e-8 idaAnalytic_mels
  TEST_ARGS ida.scalar_tolerances 1e-3 1e-8
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE idaAnalytic_mels_ida.scalar_tolerances_1e-3_1e-8.out
  EXAMPLE_TYPE exclude-single)
