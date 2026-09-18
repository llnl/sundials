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
  if(SUNDIALS_ENABLE_SUNLINSOL_KLU)
    add_subdirectory(c-klu)
  endif()
  if(SUNDIALS_ENABLE_SUNLINSOL_SUPERLUMT)
    add_subdirectory(c-superlumt)
  endif()
endif()

if(SUNDIALS_ENABLE_CXX_EXAMPLES AND SUNDIALS_ENABLE_SUNLINSOL_ONEMKLDENSE)
  add_subdirectory(cpp-onemkl)
endif()

if(SUNDIALS_ENABLE_CUDA_EXAMPLES AND SUNDIALS_ENABLE_CUDA
   AND CMAKE_CUDA_COMPILER AND SUNDIALS_INDEX_SIZE MATCHES "32")
  add_subdirectory(cuda-cusolversp)
endif()

if(SUNDIALS_ENABLE_FORTRAN_EXAMPLES AND SUNDIALS_INDEX_SIZE MATCHES "64")
  add_subdirectory(fortran)
  if(SUNDIALS_ENABLE_SUNLINSOL_LAPACKBAND
     AND SUNDIALS_ENABLE_SUNLINSOL_LAPACKDENSE)
    add_subdirectory(fortran-lapack)
  endif()
  if(SUNDIALS_ENABLE_SUNLINSOL_KLU)
    add_subdirectory(fortran-klu)
  endif()
endif()
