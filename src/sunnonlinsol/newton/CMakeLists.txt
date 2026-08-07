# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
# ------------------------------------------------------------------------------
# CMakeLists.txt file for the Newton SUNNonlinearSolver library
# ------------------------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNNONLINSOL_NEWTON\n\")")

# Add the library
sundials_add_library(
  sundials_sunnonlinsolnewton
  SOURCES sunnonlinsol_newton.c
  HEADERS ${SUNDIALS_SOURCE_DIR}/include/sunnonlinsol/sunnonlinsol_newton.h
  INCLUDE_SUBDIR sunnonlinsol
  LINK_LIBRARIES PUBLIC sundials_core
  OBJECT_LIBRARIES
  OUTPUT_NAME sundials_sunnonlinsolnewton
  VERSION ${sunnonlinsollib_VERSION}
  SOVERSION ${sunnonlinsollib_SOVERSION})

message(STATUS "Added SUNNONLINSOL_NEWTON module")

# Add F2003 module if interface is enabled
if(SUNDIALS_ENABLE_FORTRAN)
  add_subdirectory("fmod_int${SUNDIALS_INDEX_SIZE}")
endif()
