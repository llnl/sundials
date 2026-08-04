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
# CMakeLists.txt file for sunnonlinsol auto unit tests
# ------------------------------------------------------------------------------

set(examples "test_sunnonlinsol_auto\\;\\;")

include_directories(.)

set(SUNDIALS_LIBS sundials_nvecserial)
list(APPEND SUNDIALS_LIBS sundials_sunnonlinsolauto)
list(APPEND SUNDIALS_LIBS ${EXE_EXTRA_LINK_LIBS})

foreach(example_tuple ${examples})
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  if(NOT TARGET ${example})
    sundials_add_executable(${example} ${example}.c)
    set_target_properties(${example} PROPERTIES FOLDER "Examples")
    target_link_libraries(${example} ${SUNDIALS_LIBS})
  endif()

  if("${example_args}" STREQUAL "")
    set(test_name ${example})
  else()
    string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
  endif()

  sundials_add_test(
    ${test_name} ${example}
    TEST_ARGS ${example_args}
    EXAMPLE_TYPE ${example_type}
    NODIFF)
endforeach()
