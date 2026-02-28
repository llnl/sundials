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

# ---------------------------------------------------------------
# Option for testing deprecated options
# ---------------------------------------------------------------

if(DEFINED ENV{SUNDIALS_ENABLE_UNSET_DEPRECATED})
  set(_sundials_enable_unset_deprecated_default
      "$ENV{SUNDIALS_ENABLE_UNSET_DEPRECATED}")
else()
  set(_sundials_enable_unset_deprecated_default OFF)
endif()

sundials_option(
  SUNDIALS_ENABLE_UNSET_DEPRECATED BOOL "Unset deprecated CMake options"
  ${_sundials_enable_unset_deprecated_default} ADVANCED)

if(SUNDIALS_ENABLE_UNSET_DEPRECATED)
  message(WARNING "Unsetting deprecated SUNDIALS options.")
endif()

# ---------------------------------------------------------------
# Deprecated options not handled by DEPRECATED_NAMES
# ---------------------------------------------------------------

# Deprecated with SUNDIALS 6.4.0
if(DEFINED SUPERLUDIST_LIBRARY_DIR)
  message(DEPRECATION "The CMake option SUPERLUDIST_LIBRARY_DIR is deprecated. "
                      "Use SUPERLUDIST_DIR instead.")
  set(SUPERLUDIST_DIR
      "${SUPERLUDIST_LIBRARY_DIR}/../"
      CACHE BOOL "SuperLU_DIST root directory" FORCE)
  if(SUNDIALS_ENABLE_UNSET_DEPRECATED)
    unset(SUPERLUDIST_LIBRARY_DIR CACHE)
  endif()
endif()

# Deprecated CUDA_ARCH option
if(DEFINED CUDA_ARCH)
  message(DEPRECATION "The CMake option CUDA_ARCH is deprecated. "
                      "Use CMAKE_CUDA_ARCHITECTURES instead.")
  # convert sm_** to just **
  string(REGEX MATCH "[0-9]+" arch_name "${CUDA_ARCH}")
  set(CMAKE_CUDA_ARCHITECTURES
      ${arch_name}
      CACHE STRING "CUDA Architectures" FORCE)
  if(SUNDIALS_ENABLE_UNSET_DEPRECATED)
    unset(CUDA_ARCH)
  endif()
endif()

if(SUNDIALS_CALIPER_OUTPUT_DIR)
  message(
    DEPRECATION
      "The CMake option SUNDIALS_CALIPER_OUTPUT_DIR is deprecated. "
      "Use SUNDIALS_TEST_CALIPER_OUTPUT_DIR and SUNDIALS_BENCHMARK_CALIPER_OUTPUT_DIR instead."
  )
  set(SUNDIALS_TEST_CALIPER_OUTPUT_DIR
      ${SUNDIALS_TEST_CALIPER_OUTPUT_DIR}
      CACHE PATH "Location to write test caliper files" FORCE)
  set(SUNDIALS_BENCHMARK_CALIPER_OUTPUT_DIR
      ${SUNDIALS_BENCHMARK_CALIPER_OUTPUT_DIR}
      CACHE PATH "Location to write benchmark caliper files" FORCE)
  if(SUNDIALS_ENABLE_UNSET_DEPRECATED)
    unset(SUNDIALS_CALIPER_OUTPUT_DIR)
  endif()
endif()
