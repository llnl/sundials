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
# benchmarks level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

if(NOT (CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL
                                               "RelWithDebInfo"))
  message(WARNING "SUNDIALS is not being built in a \"Release\" mode, "
                  "benchmark performance will be affected")
endif()

sundials_option(BENCHMARK_NVECTOR BOOL "NVector benchmarks are on" ON)

# Disable some warnings for benchmarks
if(ENABLE_ALL_WARNINGS)
  set(CMAKE_C_FLAGS
      "${CMAKE_C_FLAGS} -Wno-unused-parameter -Wno-unused-function")
  set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -Wno-unused-parameter -Wno-unused-function")
endif()

# ----------------------------------------
# Add specific benchmarks
# ----------------------------------------

if(ENABLE_MPI)
  add_subdirectory(diffusion_2D)
  add_subdirectory(advection_reaction_3D)
endif()

# Add the nvector benchmarks
if(BENCHMARK_NVECTOR)
  add_subdirectory(nvector)
endif()
