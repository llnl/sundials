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
  idasAnalytic_mels idasAnalytic_mels.c
  EXAMPLE_TYPE exclude-single
  INSTALL_FILES idasAnalytic_mels_idas.init_step_1e-5.out)
target_link_libraries(idasAnalytic_mels PRIVATE ${IDAS_C_LIBS}
                                                sundials_nvecserial
                                                ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(
  idasAnalytic_mels_idas.init_step_1e-5 idasAnalytic_mels
  TEST_ARGS idas.init_step 1e-5
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE idasAnalytic_mels_idas.init_step_1e-5.out
  EXAMPLE_TYPE exclude-single)
