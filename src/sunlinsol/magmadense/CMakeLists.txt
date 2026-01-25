# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_MAGMADENSE\n\")")

if(SUNDIALS_MAGMA_BACKENDS MATCHES "CUDA")
  set_source_files_properties(sunlinsol_magmadense.cpp PROPERTIES LANGUAGE CUDA)
  set(_libs_needed sundials_sunmatrixmagmadense sundials_nveccuda)
elseif(SUNDIALS_MAGMA_BACKENDS MATCHES "HIP")
  set_source_files_properties(sunlinsol_magmadense.cpp PROPERTIES LANGUAGE CXX)
  set(_libs_needed sundials_sunmatrixmagmadense sundials_nvechip hip::device)
endif()

# Add the sunlinsol_magmadense library
sundials_add_library(
  sundials_sunlinsolmagmadense
  SOURCES sunlinsol_magmadense.cpp
  HEADERS ${SUNDIALS_SOURCE_DIR}/include/sunlinsol/sunlinsol_magmadense.h
  INCLUDE_SUBDIR sunlinsol
  LINK_LIBRARIES PUBLIC sundials_core
  OBJECT_LIBRARIES
  LINK_LIBRARIES PUBLIC SUNDIALS::MAGMA ${_libs_needed}
  OUTPUT_NAME sundials_sunlinsolmagmadense
  VERSION ${sunlinsollib_VERSION}
  SOVERSION ${sunlinsollib_SOVERSION})

message(STATUS "Added SUNLINSOL_MAGMADENSE module")
