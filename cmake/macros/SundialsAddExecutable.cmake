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
# CMake function that wraps the add_executable command.
#
# It adds one extra single-value argument, SCALAR_TYPE. Otherwise
# this function behaves exactly as add_executable does.
#
# ---------------------------------------------------------------

function(sundials_add_executable NAME)

  set(options)
  set(singleValueArgs SCALAR_TYPE)
  set(multiValueArgs)

  cmake_parse_arguments(arg "${options}" "${singleValueArgs}"
                        "${multiValueArgs}" ${ARGN})

  string(TOUPPER "${arg_SCALAR_TYPE}" _scalarUpper)
  if(NOT _scalarUpper OR _scalarUpper STREQUAL SUNDIALS_SCALAR_TYPE)
    add_executable(${NAME} ${arg_UNPARSED_ARGUMENTS})
  endif()

endfunction()
