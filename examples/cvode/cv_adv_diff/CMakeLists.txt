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
  if(SUNDIALS_ENABLE_SUNLINSOL_LAPACKBAND)
    add_subdirectory(c_lapack)
  endif()
  if(SUNDIALS_ENABLE_OPENMP)
    add_subdirectory(c_openmp)
  endif()
  if(SUNDIALS_ENABLE_OPENMP_DEVICE AND OPENMP_SUPPORTS_DEVICE_OFFLOADING)
    add_subdirectory(c_openmp_device)
  endif()
  if(SUNDIALS_ENABLE_MPI)
    add_subdirectory(c_mpi_diag)
    add_subdirectory(c_mpi_non)
  endif()
  if(SUNDIALS_ENABLE_HYPRE)
    add_subdirectory(c_mpi_hypre)
  endif()
  if(SUNDIALS_ENABLE_PETSC)
    add_subdirectory(c_mpi_petsc)
  endif()
endif()

if(SUNDIALS_ENABLE_CXX_EXAMPLES)
  if(SUNDIALS_ENABLE_RAJA)
    add_subdirectory(cpp_raja)
  endif()
  if(SUNDIALS_ENABLE_MPI AND SUNDIALS_ENABLE_SUPERLUDIST)
    add_subdirectory(cpp_mpi_superludist)
  endif()
endif()

if(SUNDIALS_ENABLE_CUDA_EXAMPLES)
  add_subdirectory(cuda_diag)
  add_subdirectory(cuda_kry)
  add_subdirectory(cuda_managed)
endif()

if(SUNDIALS_ENABLE_CXX_EXAMPLES AND SUNDIALS_ENABLE_HIP)
  add_subdirectory(hip_diag)
  add_subdirectory(hip_kry)
endif()

if(SUNDIALS_ENABLE_CXX_EXAMPLES AND SUNDIALS_ENABLE_SYCL)
  add_subdirectory(sycl)
endif()

if(SUNDIALS_ENABLE_FORTRAN_EXAMPLES)
  add_subdirectory(fortran)
endif()
