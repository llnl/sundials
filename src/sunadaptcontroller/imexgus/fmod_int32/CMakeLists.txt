# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ UMBC
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

# Disable warnings from unused SWIG macros which we can't control
if(ENABLE_ALL_WARNINGS)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-macros")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-macros")
endif()

sundials_add_f2003_library(
  sundials_fsunadaptcontrollerimexgus_mod
  SOURCES fsunadaptcontroller_imexgus_mod.f90 fsunadaptcontroller_imexgus_mod.c
  LINK_LIBRARIES PUBLIC sundials_fcore_mod
  OUTPUT_NAME sundials_fsunadaptcontrollerimexgus_mod OBJECT_LIB_ONLY)

message(STATUS "Added SUNAdaptController_ImExGus F2003 interface")
