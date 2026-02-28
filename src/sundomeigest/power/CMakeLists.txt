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
# CMakeLists.txt file for the PI SUNDomEigEstimator library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNDOMEIGEST_POWER\n\")")

# Add the sundomeigest_power library
sundials_add_library(
  sundials_sundomeigestpower
  SOURCES sundomeigest_power.c
  HEADERS ${SUNDIALS_SOURCE_DIR}/include/sundomeigest/sundomeigest_power.h
  INCLUDE_SUBDIR sundomeigest
  LINK_LIBRARIES PUBLIC sundials_core
  OBJECT_LIBRARIES
  OUTPUT_NAME sundials_sundomeigestpower
  VERSION ${sundomeigestlib_VERSION}
  SOVERSION ${sundomeigestlib_SOVERSION})

message(STATUS "Added SUNDOMEIGEST_POWER module")

# Add F2003 module if the interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory("fmod_int${SUNDIALS_INDEX_SIZE}")
endif()
