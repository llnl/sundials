# ---------------------------------------------------------------
# Programmer(s): Slaven Peles and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the cuda NVECTOR library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall NVECTOR_CUDA\n\")")

# Create the library
sundials_add_library(
  sundials_nveccuda
  SOURCES nvector_cuda.cu
  HEADERS ${SUNDIALS_SOURCE_DIR}/include/nvector/nvector_cuda.h
  INCLUDE_SUBDIR nvector
  LINK_LIBRARIES PUBLIC sundials_core
  OBJECT_LIBRARIES sundials_sunmemcuda_obj
  OUTPUT_NAME sundials_nveccuda
  VERSION ${nveclib_VERSION}
  SOVERSION ${nveclib_SOVERSION})

message(STATUS "Added NVECTOR_CUDA module")
