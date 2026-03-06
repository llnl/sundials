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
# CMakeLists.txt file for the sparse SUNMatrix library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNMATRIX_SPARSE\n\")")

# Add the sunmatrix_sparse library
sundials_add_library(
  sundials_sunmatrixsparse
  SOURCES sunmatrix_sparse.c
  HEADERS ${SUNDIALS_SOURCE_DIR}/include/sunmatrix/sunmatrix_sparse.h
  INCLUDE_SUBDIR sunmatrix
  LINK_LIBRARIES PUBLIC sundials_core
  OBJECT_LIBRARIES
  OUTPUT_NAME sundials_sunmatrixsparse
  VERSION ${sunmatrixlib_VERSION}
  SOVERSION ${sunmatrixlib_SOVERSION})

message(STATUS "Added SUNMATRIX_SPARSE module")

# Add F2003 module if the interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory("fmod_int${SUNDIALS_INDEX_SIZE}")
endif()
