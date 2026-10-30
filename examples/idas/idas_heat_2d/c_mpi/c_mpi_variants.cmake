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
  idas_heat_2d_fsa_kry_bbd_p SOURCES idas_heat_2d_fsa_kry_bbd_p.c
  TEST_NAME idas_heat_2d_fsa_kry_bbd_p_-sensi_stg_t
  TEST_ARGS -sensi stg t
  MPI_NPROCS 4
  EXAMPLE_TYPE exclude-single
  LINK_LIBRARIES PRIVATE MPI::MPI_C ${IDAS_C_LIBS} sundials_nvecparallel
                 ${EXE_EXTRA_LINK_LIBS})

foreach(example idas_heat_2d_kry_bbd_p idas_heat_2d_kry_p)
  sundials_add_example(
    ${example} SOURCES ${example}.c
    MPI_NPROCS 4
    EXAMPLE_TYPE exclude-single
    LINK_LIBRARIES PRIVATE MPI::MPI_C ${IDAS_C_LIBS} sundials_nvecparallel
                   ${EXE_EXTRA_LINK_LIBS})
endforeach()
