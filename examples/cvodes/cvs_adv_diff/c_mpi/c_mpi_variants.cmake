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
  cvs_adv_diff_asap_non_p
  SOURCES cvs_adv_diff_asap_non_p.c
  MPI_NPROCS 2
  EXAMPLE_TYPE exclude-single
  LINK_LIBRARIES PRIVATE MPI::MPI_C ${CVODES_C_LIBS} sundials_nvecparallel
                 ${EXE_EXTRA_LINK_LIBS})

sundials_add_example(
  cvs_adv_diff_fsa_non_p
  SOURCES cvs_adv_diff_fsa_non_p.c
  TEST_NAME cvs_adv_diff_fsa_non_p_-sensi_stg_t
  TEST_ARGS -sensi stg t
  MPI_NPROCS 2
  EXAMPLE_TYPE exclude-single
  INSTALL_FILES cvs_adv_diff_fsa_non_p_-sensi_sim_t.out
  LINK_LIBRARIES PRIVATE MPI::MPI_C ${CVODES_C_LIBS} sundials_nvecparallel
                 ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(
  cvs_adv_diff_fsa_non_p_-sensi_sim_t cvs_adv_diff_fsa_non_p
  TEST_ARGS -sensi sim t
  MPI_NPROCS 2
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE cvs_adv_diff_fsa_non_p_-sensi_sim_t.out
  EXAMPLE_TYPE exclude-single)

sundials_add_example(
  cvs_adv_diff_non_p
  SOURCES cvs_adv_diff_non_p.c
  MPI_NPROCS 2
  EXAMPLE_TYPE exclude-single
  LINK_LIBRARIES PRIVATE MPI::MPI_C ${CVODES_C_LIBS} sundials_nvecparallel
                 ${EXE_EXTRA_LINK_LIBS})
