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
# CMakeLists.txt file for the MPIPlusX NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_MPIPLUSX\n\")")

if(MPI_C_FOUND AND MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
elseif(MPI_C_FOUND)
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Create the library
sundials_add_library(
  sundials_nvecmpiplusx
  SOURCES nvector_mpiplusx.c
  HEADERS ${SUNDIALS_SOURCE_DIR}/include/nvector/nvector_mpiplusx.h
  INCLUDE_SUBDIR nvector
  LINK_LIBRARIES PUBLIC sundials_core
  OBJECT_LIBRARIES sundials_nvecmpimanyvector_obj
  OUTPUT_NAME sundials_nvecmpiplusx
  VERSION ${nveclib_VERSION}
  SOVERSION ${nveclib_SOVERSION})

message(STATUS "Added NVECTOR_MPIPLUSX module")

# Add F2003 module if the interface is enabled
if(SUNDIALS_ENABLE_FORTRAN)
  add_subdirectory("fmod_int${SUNDIALS_INDEX_SIZE}")
endif()
