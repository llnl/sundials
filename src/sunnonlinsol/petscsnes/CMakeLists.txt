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
# ------------------------------------------------------------------------------
# CMakeLists.txt file for the PETSc SUNNonlinearSolver library
# ------------------------------------------------------------------------------
install(CODE "MESSAGE(\"\nInstall SUNNONLINSOL_PETSCSNES\n\")")

if(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
else()
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Create the library
sundials_add_library(
  sundials_sunnonlinsolpetscsnes
  SOURCES sunnonlinsol_petscsnes.c
  HEADERS ${SUNDIALS_SOURCE_DIR}/include/sunnonlinsol/sunnonlinsol_petscsnes.h
  INCLUDE_SUBDIR sunnonlinsol
  LINK_LIBRARIES PUBLIC sundials_core
  OBJECT_LIBRARIES
  LINK_LIBRARIES PUBLIC sundials_nvecpetsc SUNDIALS::PETSC
  OUTPUT_NAME sundials_sunnonlinsolpetscsnes
  VERSION ${sunnonlinsollib_VERSION}
  SOVERSION ${sunnonlinsollib_SOVERSION})

message(STATUS "Added SUNNONLINSOL_PETSCSNES module")
