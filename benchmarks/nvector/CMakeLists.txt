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
# benchmarks/nvector level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

add_subdirectory(serial)

if(BUILD_NVECTOR_PARALLEL)
  add_subdirectory(parallel)
endif()

if(BUILD_NVECTOR_MPIPLUSX)
  add_subdirectory(mpiplusx)
endif()

if(BUILD_NVECTOR_PARHYP)
  add_subdirectory(parhyp)
endif()

if(BUILD_NVECTOR_PTHREADS)
  add_subdirectory(pthreads)
endif()

if(BUILD_NVECTOR_OPENMP)
  add_subdirectory(openmp)
endif()

if(BUILD_NVECTOR_OPENMPDEV)
  add_subdirectory(openmpdev)
endif()

if(BUILD_NVECTOR_PETSC)
  add_subdirectory(petsc)
endif()

if(BUILD_NVECTOR_CUDA)
  add_subdirectory(cuda)
endif()

if(BUILD_NVECTOR_HIP)
  add_subdirectory(hip)
endif()

if(BUILD_NVECTOR_RAJA)
  add_subdirectory(raja)
endif()

if(BUILD_NVECTOR_SYCL)
  add_subdirectory(sycl)
endif()

if(BUILD_NVECTOR_KOKKOS)
  add_subdirectory(kokkos)
endif()

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/plot_nvector_performance_results.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
