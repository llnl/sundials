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

# Create a library out of the generic sundials modules
sundials_add_library(
  sundials_sunmemhip
  SOURCES sundials_hip_memory.hip.cpp
  HEADERS ${SUNDIALS_SOURCE_DIR}/include/sunmemory/sunmemory_hip.h
  INCLUDE_SUBDIR sunmemory
  LINK_LIBRARIES PUBLIC sundials_core PRIVATE hip::device
  OBJECT_LIB_ONLY)
