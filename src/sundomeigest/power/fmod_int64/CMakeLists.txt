# ---------------------------------------------------------------
# Programmer(s): Mustafa Aggul @ SMU
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
# CMakeLists.txt file for the F2003 POWER SUNDOMEIGEST object library
# ---------------------------------------------------------------

sundials_add_f2003_library(
  sundials_fsundomeigestpower_mod
  SOURCES fsundomeigest_power_mod.f90 fsundomeigest_power_mod.c
  LINK_LIBRARIES PUBLIC sundials_fcore_mod
  OBJECT_LIBRARIES
  OUTPUT_NAME sundials_fsundomeigestpower_mod
  VERSION ${sundomeigestlib_VERSION}
  SOVERSION ${sundomeigestlib_SOVERSION})

message(STATUS "Added SUNDOMEIGEST_POWER F2003 interface")
