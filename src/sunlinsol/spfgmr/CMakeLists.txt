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
# CMakeLists.txt file for the SPFGMR SUNLinearSolver library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_SPFGMR\n\")")

# Add the sunlinsol_spfgmr library
sundials_add_library(
  sundials_sunlinsolspfgmr
  SOURCES sunlinsol_spfgmr.c
  HEADERS ${SUNDIALS_SOURCE_DIR}/include/sunlinsol/sunlinsol_spfgmr.h
  INCLUDE_SUBDIR sunlinsol
  LINK_LIBRARIES PUBLIC sundials_core
  OBJECT_LIBRARIES
  OUTPUT_NAME sundials_sunlinsolspfgmr
  VERSION ${sunlinsollib_VERSION}
  SOVERSION ${sunlinsollib_SOVERSION})

message(STATUS "Added SUNLINSOL_SPFGMR module")

# Add F90 module if F2003 interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory("fmod_int${SUNDIALS_INDEX_SIZE}")
endif()
