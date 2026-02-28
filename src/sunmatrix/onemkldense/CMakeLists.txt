# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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

install(CODE "MESSAGE(\"\nInstall SUNMATRIX_ONEMKLDENSE\n\")")

# Create the library
sundials_add_library(
  sundials_sunmatrixonemkldense
  SOURCES sunmatrix_onemkldense.cpp
  HEADERS ${SUNDIALS_SOURCE_DIR}/include/sunmatrix/sunmatrix_onemkldense.h
  INCLUDE_SUBDIR sunmatrix
  LINK_LIBRARIES PUBLIC sundials_core
  OBJECT_LIBRARIES sundials_sunmemsycl_obj
  COMPILE_FEATURES PUBLIC cxx_std_17
  INCLUDE_DIRECTORIES PUBLIC ${MKL_INCLUDE_DIR}
  LINK_LIBRARIES PUBLIC MKL::MKL_DPCPP
  OUTPUT_NAME sundials_sunmatrixonemkldense
  VERSION ${sunmatrixlib_VERSION}
  SOVERSION ${sunmatrixlib_SOVERSION})

message(STATUS "Added SUNMATRIX_ONEMKLDENSE module")
