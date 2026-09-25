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
  if(SUNDIALS_ENABLE_SUNLINSOL_LAPACKBAND
     AND SUNDIALS_ENABLE_SUNLINSOL_LAPACKDENSE)
    add_subdirectory(c-lapack)
  endif()
  if(SUNDIALS_ENABLE_OPENMP)
    add_subdirectory(c-openmp)
  endif()
  if(SUNDIALS_ENABLE_OPENMP_DEVICE AND OPENMP_SUPPORTS_DEVICE_OFFLOADING)
    add_subdirectory(c-openmp-device)
  endif()
  if(SUNDIALS_ENABLE_MPI)
    add_subdirectory(c-mpi)
  endif()
  if(SUNDIALS_ENABLE_HYPRE)
    add_subdirectory(c-mpi-hypre)
  endif()
  if(SUNDIALS_ENABLE_PETSC)
    add_subdirectory(c-mpi-petsc)
  endif()
endif()

if(SUNDIALS_ENABLE_CXX_EXAMPLES)
  if(SUNDIALS_ENABLE_RAJA)
    add_subdirectory(cpp-raja)
  endif()
  if(SUNDIALS_ENABLE_MPI AND SUNDIALS_ENABLE_SUPERLUDIST)
    add_subdirectory(cpp-mpi-superludist)
  endif()
endif()

if(SUNDIALS_ENABLE_CUDA_EXAMPLES AND SUNDIALS_ENABLE_CUDA
   AND CMAKE_CUDA_COMPILER)
  add_subdirectory(cuda)
  add_subdirectory(cuda-managed)
endif()

if(SUNDIALS_ENABLE_CXX_EXAMPLES AND SUNDIALS_ENABLE_HIP)
  add_subdirectory(hip)
endif()

if(SUNDIALS_ENABLE_CXX_EXAMPLES AND SUNDIALS_ENABLE_NVECTOR_SYCL)
  add_subdirectory(sycl)
endif()

if(SUNDIALS_ENABLE_FORTRAN_EXAMPLES)
  add_subdirectory(fortran)
endif()
