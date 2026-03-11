# -----------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# -----------------------------------------------------------------------------
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
# -----------------------------------------------------------------------------
# Module to find and setup oneMKL.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

include_guard(GLOBAL)

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# oneMKL does not support extended precision
if(SUNDIALS_PRECISION MATCHES "EXTENDED")
  message(
    FATAL_ERROR "oneMKL is not compatible with ${SUNDIALS_PRECISION} precision")
endif()

# oneMKL does not support 32-bit index sizes
if(SUNDIALS_INDEX_SIZE MATCHES "32")
  message(
    FATAL_ERROR
      "oneMKL is not compatible with ${SUNDIALS_INDEX_SIZE}-bit indices")
endif()

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

# Workaround bug in MKLConfig.cmake when using icpx -fsycl instead of dpcpp as
# the C++ compiler
if(SUNDIALS_ENABLE_SYCL)
  set(DPCPP_COMPILER ON)
endif()

# Look for CMake configuration file in oneMKL installation
find_package(MKL CONFIG PATHS "${ONEMKL_DIR}" "${ONEMKL_DIR}/lib/cmake/mkl"
             REQUIRED)

message(STATUS "MKL Version: ${MKL_VERSION}")
message(STATUS "MKL Targets: ${MKL_IMPORTED_TARGETS}")

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

if(SUNDIALS_ENABLE_ONEMKL_CHECKS)
  message(CHECK_START "Testing oneMKL")
  message(CHECK_PASS "success")
else()
  message(STATUS "Skipped oneMKL checks.")
endif()
