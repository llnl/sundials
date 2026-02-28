# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
# Enable SUNDIALS Benchmarks
# ---------------------------------------------------------------

#
# Check if the test runner is needed
#
if(NOT TESTRUNNER)
  # Python is needed to use the test runner
  find_package(Python3 REQUIRED)
  # Look for the testRunner script in the test directory
  find_program(
    TESTRUNNER testRunner
    PATHS test
    NO_DEFAULT_PATH REQUIRED)
  message(STATUS "Found testRunner: ${TESTRUNNER}")
  set(TESTRUNNER
      ${TESTRUNNER}
      CACHE INTERNAL "")
endif()

# Create the benchmark output directory
if(NOT EXISTS ${SUNDIALS_BENCHMARK_OUTPUT_DIR})
  file(MAKE_DIRECTORY ${SUNDIALS_BENCHMARK_OUTPUT_DIR})
endif()
message(STATUS "Benchmark output directory: ${SUNDIALS_BENCHMARK_OUTPUT_DIR}")

if(ENABLE_CALIPER)
  message(STATUS "Enabled benchmark profiling with Caliper")
  if(NOT EXISTS ${SUNDIALS_BENCHMARK_CALIPER_OUTPUT_DIR})
    file(MAKE_DIRECTORY ${SUNDIALS_BENCHMARK_CALIPER_OUTPUT_DIR})
  endif()
  message(
    STATUS
      "Benchmark Caliper output directory: ${SUNDIALS_BENCHMARK_CALIPER_OUTPUT_DIR}"
  )
endif()

#
# Create `make benchmark`
#
add_custom_target(benchmark ${CMAKE_COMMAND} -E cmake_echo_color --cyan
                            "All benchmarks complete.")
