# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# SUNMemory level CMakeLists.txt for SUNDIALS
# ------------------------------------------------------------------------------

if(SUNDIALS_ENABLE_CALIPER AND SUNDIALS_BUILD_WITH_PROFILING)
  include_directories(${caliper_INCLUDE_DIR})
endif()

add_subdirectory(system)

if(SUNDIALS_ENABLE_CUDA)
  add_subdirectory(cuda)
endif()

if(SUNDIALS_ENABLE_HIP)
  add_subdirectory(hip)
endif()

if(SUNDIALS_ENABLE_SYCL)
  add_subdirectory(sycl)
endif()
