# ---------------------------------------------------------------
# Programmer(s): Slaven Peles, Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the Trilinos NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_Trilinos\n\")")

# Create the library
sundials_add_library(
  sundials_nvectrilinos
  SOURCES nvector_trilinos.cpp
  HEADERS ${SUNDIALS_SOURCE_DIR}/include/nvector/nvector_trilinos.h
  INCLUDE_SUBDIR nvector
  LINK_LIBRARIES PUBLIC sundials_core Tpetra::all_libs
  OUTPUT_NAME sundials_nvectrilinos
  VERSION ${nveclib_VERSION}
  SOVERSION ${nveclib_SOVERSION})

# Install the exported NVECTRILINOS CXX header files in a separate subdirectory
install(
  FILES
    ${SUNDIALS_SOURCE_DIR}/include/nvector/trilinos/SundialsTpetraVectorInterface.hpp
    ${SUNDIALS_SOURCE_DIR}/include/nvector/trilinos/SundialsTpetraVectorKernels.hpp
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/nvector/trilinos")

message(STATUS "Added NVECTOR_Trilinos module")
