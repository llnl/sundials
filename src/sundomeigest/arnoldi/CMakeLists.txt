# ---------------------------------------------------------------
# Programmer(s): Mustafa Aggul @ SMU
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
# CMakeLists.txt file for the Arnoldi Iteration
# SUNDomEigEstimator library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNDOMEIGEST_ARNOLDI\n\")")

# Workaround bug in older versions of CMake where the BLAS::BLAS target, which
# LAPACK::LAPACK depends on, is not picked up by LINK_LIBRARIES
set(_lapack_targets LAPACK::LAPACK)
if(CMAKE_VERSION VERSION_LESS 3.20)
  list(APPEND _lapack_targets BLAS::BLAS)
endif()

# Add the sundomeigest_arnoldi library
sundials_add_library(
  sundials_sundomeigestarnoldi
  SOURCES sundomeigest_arnoldi.c
  HEADERS ${SUNDIALS_SOURCE_DIR}/include/sundomeigest/sundomeigest_arnoldi.h
  INCLUDE_SUBDIR sundomeigest
  LINK_LIBRARIES PUBLIC sundials_core ${_lapack_targets}
  OBJECT_LIBRARIES
  OUTPUT_NAME sundials_sundomeigestarnoldi
  VERSION ${sundomeigestlib_VERSION}
  SOVERSION ${sundomeigestlib_SOVERSION})

message(STATUS "Added SUNDOMEIGEST_ARNOLDI module")

# Add F2003 module if the interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory("fmod_int${SUNDIALS_INDEX_SIZE}")
endif()
