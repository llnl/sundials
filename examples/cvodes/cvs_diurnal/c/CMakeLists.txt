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
  cvsDiurnal_FSA_kry cvsDiurnal_FSA_kry.c
  TEST_NAME cvsDiurnal_FSA_kry_-sensi_sim_t
  TEST_ARGS -sensi sim t
  EXAMPLE_TYPE exclude-single
  INSTALL_FILES cvsDiurnal_FSA_kry_-sensi_stg_t.out
  LINK_LIBRARIES PRIVATE ${CVODES_C_LIBS} sundials_nvecserial
                 ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(
  cvsDiurnal_FSA_kry_-sensi_stg_t cvsDiurnal_FSA_kry
  TEST_ARGS -sensi stg t
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE cvsDiurnal_FSA_kry_-sensi_stg_t.out
  EXAMPLE_TYPE exclude-single)

foreach(example cvsDiurnal_kry cvsDiurnal_kry_bp)
  sundials_add_example(${example} ${example}.c)
  target_link_libraries(${example} PRIVATE ${CVODES_C_LIBS}
                                           sundials_nvecserial
                                           ${EXE_EXTRA_LINK_LIBS})
endforeach()
