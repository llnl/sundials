# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the F2003 manyvector NVECTOR object library
# ---------------------------------------------------------------

# Disable warnings from unused SWIG macros which we can't control
if(SUNDIALS_ENABLE_ALL_WARNINGS)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-macros")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-macros")
endif()

if(BUILD_NVECTOR_MANYVECTOR)
  sundials_add_f2003_library(
    sundials_fnvecmanyvector_mod
    SOURCES fnvector_manyvector_mod.f90 fnvector_manyvector_mod.c
    LINK_LIBRARIES PUBLIC sundials_fcore_mod
    OUTPUT_NAME sundials_fnvecmanyvector_mod
    VERSION ${nveclib_VERSION}
    SOVERSION ${nveclib_SOVERSION})
  message(STATUS "Added NVECTOR_MANYVECTOR F2003 Interface")
endif()

if(BUILD_NVECTOR_MPIMANYVECTOR)

  if(MPI_C_COMPILER)
    # use MPI wrapper as the compiler
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
  elseif()
    # add MPI_INCLUDE_PATH to include directories
    include_directories(${MPI_INCLUDE_PATH})
  endif()
  if(MPI_Fortran_COMPILER)
    # use MPI wrapper as the compiler
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  endif()

  sundials_add_f2003_library(
    sundials_fnvecmpimanyvector_mod
    SOURCES fnvector_mpimanyvector_mod.f90 fnvector_mpimanyvector_mod.c
    LINK_LIBRARIES PUBLIC sundials_fcore_mod
    OUTPUT_NAME sundials_fnvecmpimanyvector_mod
    VERSION ${nveclib_VERSION}
    SOVERSION ${nveclib_SOVERSION})

  message(STATUS "Added NVECTOR_MPIMANYVECTOR F2003 Interface")
endif()
