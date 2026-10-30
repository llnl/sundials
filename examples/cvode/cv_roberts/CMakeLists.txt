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
  add_subdirectory(c_dns)
  add_subdirectory(c_dns_constraints)
  add_subdirectory(c_dns_negsol)
  add_subdirectory(c_dns_uw)
  if(SUNDIALS_ENABLE_SUNLINSOL_LAPACKDENSE)
    add_subdirectory(c_lapack)
  endif()
  if(SUNDIALS_ENABLE_SUNLINSOL_KLU)
    add_subdirectory(c_klu_block)
    add_subdirectory(c_klu)
  endif()
  if(SUNDIALS_ENABLE_SUNLINSOL_SUPERLUMT)
    add_subdirectory(c_superlumt)
  endif()
endif()

if(SUNDIALS_ENABLE_CXX_EXAMPLES AND SUNDIALS_ENABLE_SUNLINSOL_ONEMKLDENSE)
  add_subdirectory(cpp_onemkl)
endif()

if(SUNDIALS_ENABLE_CUDA_EXAMPLES AND SUNDIALS_INDEX_SIZE MATCHES "32")
  add_subdirectory(cuda_cusolversp)
endif()

if(SUNDIALS_ENABLE_FORTRAN_EXAMPLES AND SUNDIALS_INDEX_SIZE MATCHES "64")
  add_subdirectory(fortran_dns)
  add_subdirectory(fortran_dns_constraints)
  if(SUNDIALS_ENABLE_SUNLINSOL_LAPACKDENSE)
    add_subdirectory(fortran_lapack)
  endif()
  if(SUNDIALS_ENABLE_SUNLINSOL_KLU)
    add_subdirectory(fortran_klu)
  endif()
endif()
