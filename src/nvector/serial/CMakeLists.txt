# ---------------------------------------------------------------
# Programmer(s): Radu Serban, Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the serial NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_SERIAL\n\")")

# Create the sundials_nvecserial library
sundials_add_library(
  sundials_nvecserial
  SOURCES nvector_serial.c
  HEADERS ${SUNDIALS_SOURCE_DIR}/include/nvector/nvector_serial.h
  INCLUDE_SUBDIR nvector
  LINK_LIBRARIES PUBLIC sundials_core
  OBJECT_LIBRARIES
  OUTPUT_NAME sundials_nvecserial
  VERSION ${nveclib_VERSION}
  SOVERSION ${nveclib_SOVERSION})

message(STATUS "Added NVECTOR_SERIAL module")

# Add F2003 module if the interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory("fmod_int${SUNDIALS_INDEX_SIZE}")
endif()
