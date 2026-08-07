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
# CMakeLists.txt file for the F2003 parallel NVECTOR object library
# ---------------------------------------------------------------

# Disable warnings from unused SWIG macros which we can't control
if(SUNDIALS_ENABLE_ALL_WARNINGS)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-macros")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-macros")
endif()

set(nvecparallel_SOURCES fnvector_parallel_mod.f90 fnvector_parallel_mod.c)

if(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
else()
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

if(MPI_Fortran_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
else()
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

sundials_add_f2003_library(
  sundials_fnvecparallel_mod
  SOURCES fnvector_parallel_mod.f90 fnvector_parallel_mod.c
  LINK_LIBRARIES PUBLIC sundials_fcore_mod
  OBJECT_LIBRARIES
  OUTPUT_NAME sundials_fnvecparallel_mod
  VERSION ${nveclib_VERSION}
  SOVERSION ${nveclib_SOVERSION})

message(STATUS "Added NVECTOR_PARALLEL F2003 interface")
