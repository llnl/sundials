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
  add_subdirectory(c)
  if(SUNDIALS_ENABLE_NVECTOR_MANYVECTOR)
    add_subdirectory(c-manyvector)
  endif()
  if(SUNDIALS_ENABLE_SUNLINSOL_KLU)
    add_subdirectory(c-klu)
  endif()
  if(SUNDIALS_ENABLE_OPENMP)
    add_subdirectory(c-openmp)
  endif()
  if(SUNDIALS_ENABLE_MPI)
    add_subdirectory(c-mpi)
  endif()
  if(SUNDIALS_ENABLE_SUNLINSOL_SUPERLUMT)
    add_subdirectory(c-superlumt)
  endif()
endif()
if(SUNDIALS_ENABLE_CXX_EXAMPLES)
  if(SUNDIALS_ENABLE_MPI AND SUNDIALS_ENABLE_RAJA
     AND SUNDIALS_PRECISION MATCHES "DOUBLE"
     AND SUNDIALS_INDEX_SIZE MATCHES "32")
    add_subdirectory(cpp-mpi)
  endif()
  if(SUNDIALS_ENABLE_SUPERLUDIST)
    add_subdirectory(cpp-mpi-superludist)
  endif()
endif()
if(SUNDIALS_ENABLE_FORTRAN_EXAMPLES AND SUNDIALS_INDEX_SIZE MATCHES "64")
  add_subdirectory(fortran)
  if(SUNDIALS_ENABLE_SUNLINSOL_KLU)
    add_subdirectory(fortran-klu)
  endif()
endif()
if(SUNDIALS_ENABLE_FORTRAN_EXAMPLES AND SUNDIALS_ENABLE_MPI)
  add_subdirectory(fortran-mpi)
endif()
