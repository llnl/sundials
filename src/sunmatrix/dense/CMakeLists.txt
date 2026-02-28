# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ UMBC
#                Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the dense SUNMatrix library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNMATRIX_DENSE\n\")")

# Add the sunmatrix_dense library
sundials_add_library(
  sundials_sunmatrixdense
  SOURCES sunmatrix_dense.c
  HEADERS ${SUNDIALS_SOURCE_DIR}/include/sunmatrix/sunmatrix_dense.h
  INCLUDE_SUBDIR sunmatrix
  LINK_LIBRARIES PUBLIC sundials_core
  OBJECT_LIBRARIES
  OUTPUT_NAME sundials_sunmatrixdense
  VERSION ${sunmatrixlib_VERSION}
  SOVERSION ${sunmatrixlib_SOVERSION})

message(STATUS "Added SUNMATRIX_DENSE module")

# Add F2003 module if the interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory("fmod_int${SUNDIALS_INDEX_SIZE}")
endif()
