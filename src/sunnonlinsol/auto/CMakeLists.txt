# ------------------------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ UMBC
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025, Lawrence Livermore National Security,
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
# CMakeLists.txt file for the fixed-point SUNNonlinearSolver library
# ------------------------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall SUNNONLINSOL_AUTO\n\")")

# Add the library
sundials_add_library(
  sundials_sunnonlinsolauto
  SOURCES sunnonlinsol_auto.c
  HEADERS ${SUNDIALS_SOURCE_DIR}/include/sunnonlinsol/sunnonlinsol_auto.h
  INCLUDE_SUBDIR sunnonlinsol
  LINK_LIBRARIES PUBLIC sundials_core sundials_sunnonlinsolnewton sundials_sunnonlinsolfixedpoint
  OBJECT_LIBRARIES
  OUTPUT_NAME sundials_sunnonlinsolauto
  VERSION ${sunnonlinsollib_VERSION}
  SOVERSION ${sunnonlinsollib_SOVERSION})

message(STATUS "Added SUNNONLINSOL_AUTO module")

# # Add F2003 module if the interface is enabled
# if(BUILD_FORTRAN_MODULE_INTERFACE)
#   add_subdirectory("fmod_int${SUNDIALS_INDEX_SIZE}")
# endif()
