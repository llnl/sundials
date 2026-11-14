# -----------------------------------------------------------------
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
# -----------------------------------------------------------------
# Preserve develop's single-precision test exclusion for these cases.
# Tests are registered separately below for their argument-specific answers.
set(_args
    "0 5"
    "1 4"
    "2 8 0 -10"
    "0 4 1"
    "0 -4"
    "1 -5"
    "2 -5 0 -10"
    "1 -3 0 -10 0"
    "0 3 0 -10 0"
    "2 4 0 -10 0"
    "0 4 0 -10 1 10 1"
    "0 4 0 -10 0 10 1")
file(GLOB _answer_files "${CMAKE_CURRENT_SOURCE_DIR}/ark_kpr_mt_*.out")
sundials_add_example(
  ark_kpr_mt
  SOURCES ark_kpr_mt.cpp
  NO_TEST
  INSTALL_FILES ${_answer_files} plot_sol.py
  INCLUDE_DIRECTORIES PRIVATE ${PROJECT_SOURCE_DIR}/examples/utilities
  LINK_LIBRARIES PRIVATE sundials_arkode sundials_nvecserial
                 ${EXE_EXTRA_LINK_LIBS})
foreach(_arg ${_args})
  string(REPLACE " " "_" _suffix "${_arg}")
  set(_type develop)
  if(_arg MATCHES "^0 -4|^1 -5|^1 -3")
    set(_type exclude-single)
  endif()
  sundials_add_test(
    ark_kpr_mt_${_suffix} ark_kpr_mt
    TEST_ARGS ${_arg}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ark_kpr_mt_${_suffix}.out
    EXAMPLE_TYPE ${_type})
endforeach()
