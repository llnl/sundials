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

# List of test tuples of the form "name\;args"
set(unit_tests "test_sunadjointcheckpointscheme_fixed\;")

# Add the build and install targets for each test
if(TARGET GTest::gtest_main AND TARGET GTest::gmock)
  foreach(test_tuple ${unit_tests})
    # parse the test tuple
    list(GET test_tuple 0 test)
    list(GET test_tuple 1 test_args)

    add_executable(${test} ${test}.cpp)
    target_link_libraries(
      ${test}
      PRIVATE sundials_adjointcheckpointscheme_fixed_obj
              sundials_sunmemsys_obj
              sundials_nvecserial
              sundials_nvecmanyvector
              sundials_core
              GTest::gtest_main
              GTest::gmock)

    gtest_discover_tests(${test})
  endforeach()
endif()
