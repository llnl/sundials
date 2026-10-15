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
  idaHeat2D_kry_tpetra idaHeat2D_kry_tpetra.cpp
  MPI_NPROCS 1
  EXAMPLE_TYPE develop
  LINK_LIBRARIES PRIVATE ${IDA_C_LIBS} sundials_nvectrilinos
                 ${EXE_EXTRA_LINK_LIBS})
if(MPI_CXX_FOUND)
  sundials_add_example(
    idaHeat2D_kry_p_tpetra idaHeat2D_kry_p_tpetra.cpp
    MPI_NPROCS 4
    EXAMPLE_TYPE develop
    LINK_LIBRARIES PRIVATE ${IDA_C_LIBS} sundials_nvectrilinos
                   ${EXE_EXTRA_LINK_LIBS})
endif()
