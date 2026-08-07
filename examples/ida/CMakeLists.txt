# ---------------------------------------------------------------
# Programmer(s): David J. Gardner, Slaven Peles, and
#                Cody Balos @ LLNL
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
# examples/ida level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# C examples
if(SUNDIALS_ENABLE_C_EXAMPLES)
  add_subdirectory(serial)
  if(SUNDIALS_ENABLE_OPENMP)
    add_subdirectory(C_openmp)
  endif()
  if(SUNDIALS_ENABLE_MPI)
    add_subdirectory(parallel)
  endif()
  if(SUNDIALS_ENABLE_PETSC)
    add_subdirectory(petsc)
  endif()
endif()

# C++ examples
if(SUNDIALS_ENABLE_CXX_EXAMPLES)
  if(SUNDIALS_ENABLE_RAJA)
    add_subdirectory(raja)
    if(SUNDIALS_ENABLE_MPI AND (SUNDIALS_RAJA_BACKENDS MATCHES "CUDA"))
      add_subdirectory(mpiraja)
    endif()
  endif()
  if(SUNDIALS_ENABLE_TRILINOS)
    add_subdirectory(trilinos)
  endif()
endif()

# Fortran examples
if(SUNDIALS_ENABLE_FORTRAN_EXAMPLES)
  add_subdirectory(F2003_serial)
  if(SUNDIALS_ENABLE_MPI)
    add_subdirectory(F2003_parallel)
  endif()
  if(SUNDIALS_ENABLE_OPENMP)
    add_subdirectory(F2003_openmp)
  endif()
endif()

# cuda examples
if(SUNDIALS_ENABLE_CUDA_EXAMPLES)

  if(SUNDIALS_ENABLE_CUDA AND CMAKE_CUDA_COMPILER)
    add_subdirectory(cuda)
    if(SUNDIALS_ENABLE_MPI)
      add_subdirectory(mpicuda)
    endif()
  endif()

endif()
