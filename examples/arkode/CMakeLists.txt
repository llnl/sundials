# ---------------------------------------------------------------
# Programmer(s): David J. Gardner and Slaven Peles @ LLNL
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
# examples/arkode level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# C examples
if(EXAMPLES_ENABLE_C)
  add_subdirectory(C_serial)
  if(BUILD_NVECTOR_MANYVECTOR)
    add_subdirectory(C_manyvector)
  endif()
  if(ENABLE_OPENMP AND OPENMP_FOUND)
    add_subdirectory(C_openmp)
  endif()
  if(ENABLE_OPENMP_DEVICE
     AND OPENMP_FOUND
     AND OPENMP_SUPPORTS_DEVICE_OFFLOADING)
    add_subdirectory(C_openmpdev)
  endif()
  if(ENABLE_MPI AND MPI_C_FOUND)
    add_subdirectory(C_parallel)
  endif()
  if(ENABLE_HYPRE AND HYPRE_FOUND)
    add_subdirectory(C_parhyp)
  endif()
  if(BUILD_SUNLINSOL_KLU)
    add_subdirectory(C_klu)
  endif()
  if(ENABLE_PETSC AND PETSC_FOUND)
    add_subdirectory(C_petsc)
  endif()
  if(BUILD_SUNLINSOL_SUPERLUMT)
    add_subdirectory(C_superlu-mt)
  endif()
endif()

# C++ examples
if(EXAMPLES_ENABLE_CXX)
  add_subdirectory(CXX_serial)
  if(ENABLE_LAPACK)
    add_subdirectory(CXX_lapack)
  endif()
  if(BUILD_NVECTOR_MANYVECTOR)
    add_subdirectory(CXX_manyvector)
  endif()
  if(ENABLE_MPI AND MPI_CXX_FOUND)
    add_subdirectory(CXX_parallel)
  endif()
  if(ENABLE_MPI
     AND MPI_CXX_FOUND
     AND ENABLE_HYPRE
     AND HYPRE_FOUND)
    add_subdirectory(CXX_parhyp)
  endif()
  if(MPI_CXX_FOUND AND ENABLE_SUPERLUDIST)
    add_subdirectory(CXX_superludist)
  endif()
  if(ENABLE_XBRAID)
    add_subdirectory(CXX_xbraid)
  endif()
endif()

# F2003 examples
if(BUILD_FORTRAN_MODULE_INTERFACE AND EXAMPLES_ENABLE_F2003)
  add_subdirectory(F2003_serial)
  add_subdirectory(F2003_custom)
  if(ENABLE_MPI AND MPI_Fortran_FOUND)
    add_subdirectory(F2003_parallel)
  endif()
endif()
