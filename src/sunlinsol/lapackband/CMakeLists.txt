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
# CMakeLists.txt file for the LapackBand SUNLinearSolver library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_LAPACKBAND\n\")")

# Workaround bug in older versions of CMake where the BLAS::BLAS target, which
# LAPACK::LAPACK depends on, is not picked up by LINK_LIBRARIES
set(_lapack_targets LAPACK::LAPACK)
if(CMAKE_VERSION VERSION_LESS 3.20)
  list(APPEND _lapack_targets BLAS::BLAS)
endif()

# Add the library
sundials_add_library(
  sundials_sunlinsollapackband
  SOURCES sunlinsol_lapackband.c
  HEADERS ${SUNDIALS_SOURCE_DIR}/include/sunlinsol/sunlinsol_lapackband.h
  INCLUDE_SUBDIR sunlinsol
  LINK_LIBRARIES PUBLIC sundials_core
  OBJECT_LIBRARIES
  LINK_LIBRARIES PUBLIC sundials_sunmatrixband ${_lapack_targets}
  OUTPUT_NAME sundials_sunlinsollapackband
  VERSION ${sunlinsollib_VERSION}
  SOVERSION ${sunlinsollib_SOVERSION})

message(STATUS "Added SUNLINSOL_LAPACKBAND module")
