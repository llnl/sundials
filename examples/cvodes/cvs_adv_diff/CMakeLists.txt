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
if(SUNDIALS_ENABLE_C_EXAMPLES)
  add_subdirectory(c_adv_diff_asai_bnd)
  add_subdirectory(c_adv_diff_bnd)
  add_subdirectory(c_adv_diff_fsa_non)
  if(SUNDIALS_ENABLE_SUNLINSOL_LAPACKBAND
     AND SUNDIALS_ENABLE_SUNLINSOL_LAPACKDENSE)
    add_subdirectory(c_lapack)
  endif()
  if(SUNDIALS_ENABLE_OPENMP)
    add_subdirectory(c_openmp)
  endif()
  if(SUNDIALS_ENABLE_MPI)
    add_subdirectory(c_mpi_cvs_adv_diff_asap_non_p)
    add_subdirectory(c_mpi_cvs_adv_diff_fsa_non_p)
    add_subdirectory(c_mpi_cvs_adv_diff_non_p)
  endif()
endif()

if(SUNDIALS_ENABLE_FORTRAN_EXAMPLES)
  add_subdirectory(fortran)
endif()
