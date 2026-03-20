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

install(CODE "MESSAGE(\"\nInstall SUNLINSOL_ONEMKLDENSE\n\")")

# Add the sunlinsol_onemkldense library
sundials_add_library(
  sundials_sunlinsolonemkldense
  SOURCES sunlinsol_onemkldense.cpp
  HEADERS ${SUNDIALS_SOURCE_DIR}/include/sunlinsol/sunlinsol_onemkldense.h
  INCLUDE_SUBDIR sunlinsol
  LINK_LIBRARIES PUBLIC sundials_core
  OBJECT_LIBRARIES
  COMPILE_FEATURES PUBLIC cxx_std_17
  INCLUDE_DIRECTORIES PUBLIC ${MKL_INCLUDE_DIR}
  LINK_LIBRARIES PUBLIC MKL::MKL_DPCPP sundials_sunmatrixonemkldense
                 sundials_nvecsycl
  OUTPUT_NAME sundials_sunlinsolonemkldense
  VERSION ${sunlinsollib_VERSION}
  SOVERSION ${sunlinsollib_SOVERSION})

message(STATUS "Added SUNLINSOL_ONEMKLDENSE module")
