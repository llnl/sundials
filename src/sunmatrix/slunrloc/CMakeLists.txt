# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ------------------------------------------------------------------------------
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
# -----------------------------------------------------------------------------
# CMakeLists.txt file for the SuperLU SLU_NR_LOC SUNMatrix library
# -----------------------------------------------------------------------------

install(CODE "message(\"\nInstall SUNMATRIX_SLUNRLOC\n\")")

# Use MPI
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

# Add the library
sundials_add_library(
  sundials_sunmatrixslunrloc
  SOURCES sunmatrix_slunrloc.c
  HEADERS ${SUNDIALS_SOURCE_DIR}/include/sunmatrix/sunmatrix_slunrloc.h
  INCLUDE_SUBDIR sunmatrix
  LINK_LIBRARIES PUBLIC sundials_core
  OBJECT_LIBRARIES
  LINK_LIBRARIES PUBLIC SUNDIALS::SUPERLUDIST
                 $<IF:$<BOOL:${SUPERLUDIST_OpenMP}>,OpenMP::OpenMP_C,>
  COMPILE_OPTIONS PRIVATE ${_compile_options}
  OUTPUT_NAME sundials_sunmatrixslunrloc
  VERSION ${sunlinsollib_VERSION}
  SOVERSION ${sunlinsollib_VERSION})

message(STATUS "Added SUNMATRIX_SLUNRLOC module")
