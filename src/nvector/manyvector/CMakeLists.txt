# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ UMBC,
#                and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the ManyVector NVECTOR library
# ---------------------------------------------------------------

# Create the sundials_nvecmanyvector library
if(SUNDIALS_ENABLE_NVECTOR_MANYVECTOR)
  install(CODE "MESSAGE(\"\nInstall NVECTOR_MANYVECTOR\n\")")

  sundials_add_library(
    sundials_nvecmanyvector
    SOURCES nvector_manyvector.c
    HEADERS ${SUNDIALS_SOURCE_DIR}/include/nvector/nvector_manyvector.h
    INCLUDE_SUBDIR nvector
    LINK_LIBRARIES PUBLIC sundials_core
    OUTPUT_NAME sundials_nvecmanyvector
    VERSION ${nveclib_VERSION}
    SOVERSION ${nveclib_SOVERSION})

  message(STATUS "Added NVECTOR_MANYVECTOR module")
endif()

# Create the sundials_nvecmpimanyvector library
if(SUNDIALS_ENABLE_NVECTOR_MPIMANYVECTOR)
  install(CODE "MESSAGE(\"\nInstall NVECTOR_MPIMANYVECTOR\n\")")

  if(MPI_C_COMPILER)
    # use MPI wrapper as the compiler
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
  else()
    # add MPI_INCLUDE_PATH to include directories
    include_directories(${MPI_INCLUDE_PATH})
  endif()

  # Create the sundials_nvecmpimanyvector library
  sundials_add_library(
    sundials_nvecmpimanyvector
    SOURCES nvector_manyvector.c
    HEADERS ${SUNDIALS_SOURCE_DIR}/include/nvector/nvector_mpimanyvector.h
    INCLUDE_SUBDIR nvector
    LINK_LIBRARIES PUBLIC sundials_core
    COMPILE_DEFINITIONS PRIVATE MANYVECTOR_BUILD_WITH_MPI
    OUTPUT_NAME sundials_nvecmpimanyvector
    VERSION ${nveclib_VERSION}
    SOVERSION ${nveclib_SOVERSION})

  message(STATUS "Added NVECTOR_MPIMANYVECTOR module")
endif()

# Add F2003 module if the interface is enabled
if(SUNDIALS_ENABLE_FORTRAN)
  add_subdirectory("fmod_int${SUNDIALS_INDEX_SIZE}")
endif()
