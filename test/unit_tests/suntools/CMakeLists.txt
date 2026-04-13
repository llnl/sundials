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

# Locate suntools package source directory.
set(SUNTOOLS_DIR "${CMAKE_SOURCE_DIR}/suntools/src")

if(NOT EXISTS "${SUNTOOLS_DIR}/suntools/logs.py")
  message(WARNING "Skipped suntools log parser tests (suntools/logs.py not found)")
  return()
endif()

# TODO(DJG): Update sundials_add_test to better handle python tests

# Add the test and configure the environment
add_test(NAME test_logfile_parse_payload
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_logfile_parse_payload.py)

set_tests_properties(
  test_logfile_parse_payload
  PROPERTIES ENVIRONMENT "PYTHONPATH=${SUNTOOLS_DIR}" LABELS "suntools;logging")

# Add the test and configure the environment
add_test(NAME test_logfile_parse_line
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_logfile_parse_line.py)

set_tests_properties(
  test_logfile_parse_line PROPERTIES ENVIRONMENT "PYTHONPATH=${SUNTOOLS_DIR}"
                                     LABELS "suntools;logging")

# Add the test and configure the environment
add_test(NAME test_logfile_to_list
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_logfile_to_list.py)

set_tests_properties(
  test_logfile_to_list PROPERTIES ENVIRONMENT "PYTHONPATH=${SUNTOOLS_DIR}"
                                  LABELS "suntools;logging")

# Add the test and configure the environment
add_test(NAME test_cli_parse_logs
         COMMAND ${Python3_EXECUTABLE}
                 ${CMAKE_CURRENT_SOURCE_DIR}/test_cli_parse_logs.py)

set_tests_properties(
  test_cli_parse_logs PROPERTIES ENVIRONMENT "PYTHONPATH=${SUNTOOLS_DIR}"
                                 LABELS "suntools;logging")

message(STATUS "Added suntools log parser tests")
