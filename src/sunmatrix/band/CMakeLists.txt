# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ UMBC
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
# CMakeLists.txt file for the band SUNMatrix library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNMATRIX_BAND\n\")")

# Add the sunmatrix_band library
sundials_add_library(
  sundials_sunmatrixband
  SOURCES sunmatrix_band.c
  HEADERS ${SUNDIALS_SOURCE_DIR}/include/sunmatrix/sunmatrix_band.h
  INCLUDE_SUBDIR sunmatrix
  LINK_LIBRARIES PUBLIC sundials_core
  OBJECT_LIBRARIES
  OUTPUT_NAME sundials_sunmatrixband
  VERSION ${sunmatrixlib_VERSION}
  SOVERSION ${sunmatrixlib_SOVERSION})

message(STATUS "Added SUNMATRIX_BAND module")

# Add F2003 module if the interface is enabled
if(SUNDIALS_ENABLE_FORTRAN)
  add_subdirectory("fmod_int${SUNDIALS_INDEX_SIZE}")
endif()
