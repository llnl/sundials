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
  if(SUNDIALS_ENABLE_SUNLINSOL_KLU)
    add_subdirectory(c_klu)
  endif()
  if(SUNDIALS_ENABLE_MPI)
    add_subdirectory(c_mpi)
  endif()
  if(SUNDIALS_ENABLE_PETSC)
    add_subdirectory(c_mpi_petsc)
  endif()
endif()

if(SUNDIALS_ENABLE_CXX_EXAMPLES)
  if(SUNDIALS_ENABLE_RAJA)
    add_subdirectory(cpp_raja)
    if(SUNDIALS_ENABLE_MPI AND SUNDIALS_RAJA_BACKENDS MATCHES "CUDA")
      add_subdirectory(cpp_mpi_raja)
    endif()
  endif()
  if(SUNDIALS_ENABLE_TRILINOS)
    add_subdirectory(cpp_trilinos)
  endif()
endif()

if(SUNDIALS_ENABLE_FORTRAN_EXAMPLES AND SUNDIALS_INDEX_SIZE MATCHES "64")
  add_subdirectory(fortran)
  if(SUNDIALS_ENABLE_MPI)
    add_subdirectory(fortran_mpi)
  endif()
  if(SUNDIALS_ENABLE_OPENMP)
    add_subdirectory(fortran_openmp)
  endif()
endif()

if(SUNDIALS_ENABLE_CUDA_EXAMPLES
   AND SUNDIALS_ENABLE_CUDA
   AND CMAKE_CUDA_COMPILER)
  add_subdirectory(cuda)
  if(SUNDIALS_ENABLE_MPI)
    add_subdirectory(cuda_mpi)
  endif()
endif()
