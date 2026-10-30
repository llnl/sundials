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
  cvs_diurnal_fsa_kry SOURCES cvs_diurnal_fsa_kry.c
  TEST_NAME cvs_diurnal_fsa_kry_-sensi_sim_t
  TEST_ARGS -sensi sim t
  EXAMPLE_TYPE exclude-single
  INSTALL_FILES cvs_diurnal_fsa_kry_-sensi_stg_t.out
  LINK_LIBRARIES PRIVATE ${CVODES_C_LIBS} sundials_nvecserial
                 ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(
  cvs_diurnal_fsa_kry_-sensi_stg_t cvs_diurnal_fsa_kry
  TEST_ARGS -sensi stg t
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE cvs_diurnal_fsa_kry_-sensi_stg_t.out
  EXAMPLE_TYPE exclude-single)

foreach(example cvs_diurnal_kry cvs_diurnal_kry_bp)
  sundials_add_example(
    ${example} SOURCES ${example}.c
    LINK_LIBRARIES PRIVATE ${CVODES_C_LIBS} sundials_nvecserial
                   ${EXE_EXTRA_LINK_LIBS})
endforeach()
