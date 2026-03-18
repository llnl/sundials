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
if(SUNDIALS_ENABLE_C_EXAMPLES)
  add_subdirectory(C_serial)
  if(SUNDIALS_ENABLE_NVECTOR_MANYVECTOR)
    add_subdirectory(C_manyvector)
  endif()
  if(SUNDIALS_ENABLE_OPENMP)
    add_subdirectory(C_openmp)
  endif()
  if(SUNDIALS_ENABLE_OPENMP_DEVICE AND OPENMP_SUPPORTS_DEVICE_OFFLOADING)
    add_subdirectory(C_openmpdev)
  endif()
  if(SUNDIALS_ENABLE_MPI)
    add_subdirectory(C_parallel)
  endif()
  if(SUNDIALS_ENABLE_HYPRE)
    add_subdirectory(C_parhyp)
  endif()
  if(SUNDIALS_ENABLE_SUNLINSOL_KLU)
    add_subdirectory(C_klu)
  endif()
  if(SUNDIALS_ENABLE_PETSC)
    add_subdirectory(C_petsc)
  endif()
  if(SUNDIALS_ENABLE_SUNLINSOL_SUPERLUMT)
    add_subdirectory(C_superlu-mt)
  endif()
endif()

# C++ examples
if(SUNDIALS_ENABLE_CXX_EXAMPLES)
  add_subdirectory(CXX_serial)
  if(SUNDIALS_ENABLE_LAPACK)
    add_subdirectory(CXX_lapack)
  endif()
  if(SUNDIALS_ENABLE_NVECTOR_MANYVECTOR)
    add_subdirectory(CXX_manyvector)
  endif()
  if(SUNDIALS_ENABLE_MPI)
    add_subdirectory(CXX_parallel)
  endif()
  if(SUNDIALS_ENABLE_HYPRE)
    add_subdirectory(CXX_parhyp)
  endif()
  if(SUNDIALS_ENABLE_SUPERLUDIST)
    add_subdirectory(CXX_superludist)
  endif()
  if(SUNDIALS_ENABLE_XBRAID)
    add_subdirectory(CXX_xbraid)
  endif()
endif()

# F2003 examples
if(SUNDIALS_ENABLE_FORTRAN_EXAMPLES)
  add_subdirectory(F2003_serial)
  add_subdirectory(F2003_custom)
  if(SUNDIALS_ENABLE_MPI)
    add_subdirectory(F2003_parallel)
  endif()
endif()
