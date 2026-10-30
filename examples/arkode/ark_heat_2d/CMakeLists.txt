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
if(SUNDIALS_ENABLE_CXX_EXAMPLES)
  add_subdirectory(cpp_heat_2d)
  add_subdirectory(cpp_heat_2d_lsrk)
  if(SUNDIALS_ENABLE_LAPACK)
    add_subdirectory(cpp_lapack)
  endif()
  if(SUNDIALS_ENABLE_MPI)
    add_subdirectory(cpp_mpi)
    if(SUNDIALS_ENABLE_HYPRE)
      add_subdirectory(cpp_mpi_hypre)
    endif()
    if(SUNDIALS_ENABLE_XBRAID)
      add_subdirectory(cpp_mpi_xbraid)
    endif()
  endif()
endif()
if(SUNDIALS_ENABLE_FORTRAN_EXAMPLES
   AND SUNDIALS_ENABLE_MPI
   AND SUNDIALS_INDEX_SIZE MATCHES "64")
  add_subdirectory(fortran_mpi)
endif()
