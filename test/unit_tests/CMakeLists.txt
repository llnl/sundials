# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025, Lawrence Livermore National Security,
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
# unit_tests level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

message("Adding units tests")

# Disable some warnings for tests
if(ENABLE_ALL_WARNINGS)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unused-parameter")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter")
endif()

add_subdirectory(sundials)

if(BUILD_ARKODE)
  add_subdirectory(arkode)
endif()

if(BUILD_CVODE)
  add_subdirectory(cvode)
endif()

if(BUILD_CVODES)
  add_subdirectory(cvodes)
endif()

if(BUILD_IDA)
  add_subdirectory(ida)
endif()

if(BUILD_IDAS)
  add_subdirectory(idas)
endif()

if(BUILD_KINSOL)
  add_subdirectory(kinsol)
endif()

add_subdirectory(constraints)

add_subdirectory(nvector)

add_subdirectory(sunmatrix)

add_subdirectory(sunlinsol)

add_subdirectory(sunnonlinsol)

add_subdirectory(sunadjointcheckpointscheme)

add_subdirectory(sundomeigest)

if(CXX_FOUND)
  add_subdirectory(logging)
  add_subdirectory(sunmemory)
  if(SUNDIALS_BUILD_WITH_PROFILING)
    add_subdirectory(profiling)
  endif()
endif()
