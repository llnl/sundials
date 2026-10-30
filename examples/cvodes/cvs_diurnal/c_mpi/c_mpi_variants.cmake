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
  cvs_diurnal_fsa_kry_p
  SOURCES cvs_diurnal_fsa_kry_p.c
  TEST_NAME cvs_diurnal_fsa_kry_p_-sensi_stg_t
  TEST_ARGS -sensi stg t
  MPI_NPROCS 4
  EXAMPLE_TYPE exclude-single
  INSTALL_FILES cvs_diurnal_fsa_kry_p_-sensi_sim_t.out
  LINK_LIBRARIES PRIVATE MPI::MPI_C ${CVODES_C_LIBS} sundials_nvecparallel
                 ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(
  cvs_diurnal_fsa_kry_p_-sensi_sim_t cvs_diurnal_fsa_kry_p
  TEST_ARGS -sensi sim t
  MPI_NPROCS 4
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE cvs_diurnal_fsa_kry_p_-sensi_sim_t.out
  EXAMPLE_TYPE exclude-single)

foreach(example cvs_diurnal_kry_bbd_p cvs_diurnal_kry_p)
  sundials_add_example(
    ${example}
    SOURCES ${example}.c
    MPI_NPROCS 4
    EXAMPLE_TYPE exclude-single
    LINK_LIBRARIES PRIVATE MPI::MPI_C ${CVODES_C_LIBS} sundials_nvecparallel
                   ${EXE_EXTRA_LINK_LIBS})
endforeach()
