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
  add_subdirectory(c_brusselator)
  add_subdirectory(c_brusselator_1d)
  add_subdirectory(c_brusselator_1d_imexmri)
  add_subdirectory(c_brusselator_1d_mri)
  add_subdirectory(c_brusselator_fp)
  add_subdirectory(c_brusselator_lsrk_domeigest)
  add_subdirectory(c_brusselator_lsrk_externaldomeigest)
  add_subdirectory(c_brusselator_mri)
  if(SUNDIALS_ENABLE_NVECTOR_MANYVECTOR)
    add_subdirectory(c_manyvector)
  endif()
  if(SUNDIALS_ENABLE_SUNLINSOL_KLU)
    add_subdirectory(c_klu)
  endif()
  if(SUNDIALS_ENABLE_OPENMP)
    add_subdirectory(c_openmp)
  endif()
  if(SUNDIALS_ENABLE_MPI
     AND SUNDIALS_PRECISION MATCHES "DOUBLE"
     AND SUNDIALS_INDEX_SIZE MATCHES "64")
    add_subdirectory(c_mpi)
  endif()
  if(SUNDIALS_ENABLE_SUNLINSOL_SUPERLUMT)
    add_subdirectory(c_superlumt)
  endif()
endif()
if(SUNDIALS_ENABLE_CXX_EXAMPLES)
  if(SUNDIALS_ENABLE_MPI
     AND SUNDIALS_ENABLE_RAJA
     AND SUNDIALS_PRECISION MATCHES "DOUBLE"
     AND SUNDIALS_INDEX_SIZE MATCHES "32")
    add_subdirectory(cpp_mpi_brusselator_1d_task_local_nls_use_cuda_nvec)
    add_subdirectory(cpp_mpi_brusselator_1d_task_local_nls_use_cudauvm_nvec)
    add_subdirectory(cpp_mpi_brusselator_1d_task_local_nls_use_hip_nvec)
    add_subdirectory(cpp_mpi_brusselator_1d_task_local_nls_use_ompdev_nvec)
    add_subdirectory(cpp_mpi_brusselator_1d_task_local_nls_use_raja_nvec)
    add_subdirectory(cpp_mpi_brusselator_1d_task_local_nls_use_serial_nvec)
  endif()
  if(SUNDIALS_ENABLE_SUPERLUDIST)
    add_subdirectory(cpp_mpi_superludist)
  endif()
endif()
if(SUNDIALS_ENABLE_FORTRAN_EXAMPLES AND SUNDIALS_INDEX_SIZE MATCHES "64")
  add_subdirectory(fortran)
  if(SUNDIALS_ENABLE_SUNLINSOL_KLU)
    add_subdirectory(fortran_klu)
  endif()
endif()
if(SUNDIALS_ENABLE_FORTRAN_EXAMPLES AND SUNDIALS_ENABLE_MPI)
  add_subdirectory(fortran_mpi)
endif()
