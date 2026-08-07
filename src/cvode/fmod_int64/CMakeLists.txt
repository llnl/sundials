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
# CMakeLists.txt file for the F2003 CVODE object library
# ---------------------------------------------------------------

# Disable warnings from unused SWIG macros which we can't control
if(SUNDIALS_ENABLE_ALL_WARNINGS)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-macros")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-macros")
endif()

set(cvode_SOURCES fcvode_mod.f90 fcvode_mod.c)

# Create the library
sundials_add_f2003_library(
  sundials_fcvode_mod
  SOURCES ${cvode_SOURCES}
  LINK_LIBRARIES PUBLIC sundials_fcore_mod
  OBJECT_LIBRARIES
    sundials_fnvecserial_mod_obj
    sundials_fsunmatrixband_mod_obj
    sundials_fsunmatrixdense_mod_obj
    sundials_fsunmatrixsparse_mod_obj
    sundials_fsunlinsolband_mod_obj
    sundials_fsunlinsoldense_mod_obj
    sundials_fsunlinsolspbcgs_mod_obj
    sundials_fsunlinsolspfgmr_mod_obj
    sundials_fsunlinsolspgmr_mod_obj
    sundials_fsunlinsolsptfqmr_mod_obj
    sundials_fsunlinsolpcg_mod_obj
    sundials_fsunnonlinsolnewton_mod_obj
    sundials_fsunnonlinsolfixedpoint_mod_obj
  OUTPUT_NAME sundials_fcvode_mod
  VERSION ${cvodelib_VERSION}
  SOVERSION ${cvodelib_SOVERSION})

message(STATUS "Added CVODE F2003 interface")
