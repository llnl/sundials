# ---------------------------------------------------------------
# Programmer(s): David J. Gardner, Slaven Peles, and
#                Shelby Lockhart @ LLNL
# ---------------------------------------------------------------
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
# ---------------------------------------------------------------
# examples/cvode level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# C examples
if(EXAMPLES_ENABLE_C)
  add_subdirectory(serial)
  if(SUNDIALS_ENABLE_OPENMP)
    add_subdirectory(C_openmp)
  endif()
  if(SUNDIALS_ENABLE_OPENMP_DEVICE AND OPENMP_SUPPORTS_DEVICE_OFFLOADING)
    add_subdirectory(C_openmpdev)
  endif()
  if(SUNDIALS_ENABLE_MPI)
    add_subdirectory(parallel)
    if(BUILD_NVECTOR_MPIMANYVECTOR)
      add_subdirectory(C_mpimanyvector)
    endif()
  endif()
  if(SUNDIALS_ENABLE_HYPRE)
    add_subdirectory(parhyp)
  endif()
  if(SUNDIALS_ENABLE_PETSC)
    add_subdirectory(petsc)
  endif()
endif()

# CXX examples
if(EXAMPLES_ENABLE_CXX)
  add_subdirectory(CXX_serial)
  if(SUNDIALS_ENABLE_HIP)
    add_subdirectory(hip)
  endif()
  if(SUNDIALS_ENABLE_MAGMA)
    add_subdirectory(magma)
  endif()
  if(SUNDIALS_ENABLE_GINKGO)
    add_subdirectory(ginkgo)
  endif()
  if(BUILD_SUNLINSOL_KOKKOSDENSE)
    add_subdirectory(kokkos)
  endif()
  if(BUILD_SUNLINSOL_ONEMKLDENSE)
    add_subdirectory(CXX_onemkl)
  endif()
  if(SUNDIALS_ENABLE_RAJA)
    add_subdirectory(raja)
  endif()
  if(BUILD_NVECTOR_SYCL)
    add_subdirectory(CXX_sycl)
  endif()
  # MPI based examples
  if(SUNDIALS_ENABLE_MPI)
    add_subdirectory(CXX_parallel)
    # hypre examples
    if(SUNDIALS_ENABLE_HYPRE)
      add_subdirectory(CXX_parhyp)
    endif()
    # SuperLU_DIST examples
    if(SUNDIALS_ENABLE_SUPERLUDIST)
      add_subdirectory(superludist)
    endif()
  endif()
endif()

# Fortran examples
if(EXAMPLES_ENABLE_F2003)
  add_subdirectory(F2003_serial)
  if(SUNDIALS_ENABLE_MPI)
    add_subdirectory(F2003_parallel)
  endif()
endif()

# CUDA examples
if(EXAMPLES_ENABLE_CUDA)

  if(SUNDIALS_ENABLE_CUDA AND CMAKE_CUDA_COMPILER)
    add_subdirectory(cuda)
  endif()

endif()
