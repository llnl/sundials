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

### TESTS THAT REQUIRE GTEST ###

# List of test tuples of the form "name\;args"
set(unit_tests "test_sundials_datanode\;" "test_sundials_stlvector\;"
               "test_sundials_hashmap\;" "test_sundials_logger\;")

if(SUNDIALS_ENABLE_ERROR_CHECKS)
  list(APPEND gtest_tests "test_sundials_errors\;")
endif()

# Add the build and install targets for each test
if(TARGET GTest::gtest_main AND TARGET GTest::gmock)
  foreach(test_tuple ${gtest_tests})
    # parse the test tuple
    list(GET test_tuple 0 test)
    list(GET test_tuple 1 test_args)

    add_executable(
      ${test} ${test}.cpp ${PROJECT_SOURCE_DIR}/src/sundials/sundials_hashmap.c)
    target_link_libraries(
      ${test} PRIVATE sundials_core sundials_nvecserial sundials_sunmemsys_obj
                      GTest::gtest_main GTest::gmock)

    if(test STREQUAL "test_sundials_logger")
      gtest_discover_tests(${test} PROPERTIES LABELS "logging")
    else()
      gtest_discover_tests(${test})
    endif()
  endforeach()
endif()


### REMAINING TESTS ###

# List of test tuples of the form "name\;args"
set(unit_tests "test_sundials_QRfactorization\;" "test_sundials_GramSchmidt\;"
               "test_sundials_QRAddVariants\;")

# Add the build and install targets for each test
foreach(test_tuple ${unit_tests})
  # parse the test tuple
  list(GET test_tuple 0 test)
  list(GET test_tuple 1 test_args)

  add_executable(
    ${test} ${test}.cpp)
  target_link_libraries(
    ${test} PRIVATE sundials_core sundials_nvecserial)
endforeach()
