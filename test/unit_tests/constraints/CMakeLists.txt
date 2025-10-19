# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025, Lawrence Livermore National Security,
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

# List of test tuples of the form "name"
set(unit_tests)

if(BUILD_ARKODE)
  list(APPEND unit_tests "test_constraints_arkode.cpp")
endif()

if(BUILD_CVODE)
  list(APPEND unit_tests "test_constraints_cvode.cpp")
endif()

if(BUILD_CVODES)
  list(APPEND unit_tests "test_constraints_cvodes.cpp")
endif()

if(BUILD_IDA)
  list(APPEND unit_tests "test_constraints_ida.cpp")
endif()

if(BUILD_IDAS)
  list(APPEND unit_tests "test_constraints_idas.cpp")
endif()

# Add the build and install targets for each test
foreach(test_tuple ${unit_tests})

  # parse the test tuple
  list(GET test_tuple 0 test_src)

  # extract the file name without extension
  get_filename_component(test_target ${test_src} NAME_WE)

  # The CMake regular expression support is somewhat limited (see this issue
  # https://gitlab.kitware.com/cmake/cmake/-/issues/17686) so the ordering here
  # matters i.e., match cvodes before cvode and idas before ida
  string(REGEX MATCH "arkode|cvodes|cvode|idas|ida|kinsol" _pkg ${test_target})

  # check if this test has already been added, only need to add test source
  # files once for testing with different inputs
  if(NOT TARGET ${test_target})

    # test source files
    sundials_add_executable(${test_target} ${test_src})

    set_target_properties(${test_target} PROPERTIES FOLDER "unit_tests")

    # include location of public and private header files
    target_include_directories(
      ${test_target} PRIVATE ${CMAKE_SOURCE_DIR}/include
                             ${CMAKE_SOURCE_DIR}/test/unit_tests)

    # libraries to link against
    target_link_libraries(${test_target} sundials_${_pkg} sundials_nvecserial
                          ${EXE_EXTRA_LINK_LIBS})

  endif()

  # add test to regression tests
  sundials_add_test(
    ${test_target} ${test_target}
    LABELS "constraints"
    NODIFF)

endforeach()

message(STATUS "Added constraints units tests")
