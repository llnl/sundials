# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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

# List of test tuples of the form "name\;args"
set(unit_tests)

if(BUILD_ARKODE)
  # ARKStep
  list(APPEND unit_tests "test_logging_arkode_arkstep.cpp\;0") # ERK
  list(APPEND unit_tests "test_logging_arkode_arkstep.cpp\;1 1 1") # DIRK Newton
                                                                   # + Dense
  list(APPEND unit_tests "test_logging_arkode_arkstep.cpp\;1 1 0") # DIRK Newton
                                                                   # + GMRES
  list(APPEND unit_tests "test_logging_arkode_arkstep.cpp\;1 0") # DIRK
                                                                 # Fixed-point
  list(APPEND unit_tests "test_logging_arkode_arkstep.cpp\;2 1 1") # ImEx Newton
                                                                   # + Dense
  list(APPEND unit_tests "test_logging_arkode_arkstep.cpp\;2 1 0") # ImEx Newton
                                                                   # + GMRES
  list(APPEND unit_tests "test_logging_arkode_arkstep.cpp\;2 0") # ImEx
                                                                 # Fixed-point
  # ERKStep
  list(APPEND unit_tests "test_logging_arkode_erkstep.cpp\;")
  # ForcingStep
  list(APPEND unit_tests "test_logging_arkode_forcingstep.cpp\;")
  # LSRKStep
  list(APPEND unit_tests "test_logging_arkode_lsrkstep.cpp\;0") # RKC
  list(APPEND unit_tests "test_logging_arkode_lsrkstep.cpp\;1") # RKL
  list(APPEND unit_tests "test_logging_arkode_lsrkstep.cpp\;2") # SSPs2
  list(APPEND unit_tests "test_logging_arkode_lsrkstep.cpp\;3") # SSPs3
  list(APPEND unit_tests "test_logging_arkode_lsrkstep.cpp\;4") # SSP43
  list(APPEND unit_tests "test_logging_arkode_lsrkstep.cpp\;5") # SSP104
  # MRIStep -- MRI-GARK
  list(APPEND unit_tests "test_logging_arkode_mristep.cpp\;0") # Explicit
  list(APPEND unit_tests "test_logging_arkode_mristep.cpp\;1 1 1") # Implicit
                                                                   # Newton +
                                                                   # Dense
  list(APPEND unit_tests "test_logging_arkode_mristep.cpp\;1 1 0") # Implicit
                                                                   # Newton +
                                                                   # GMRES
  list(APPEND unit_tests "test_logging_arkode_mristep.cpp\;1 0") # Implicit
                                                                 # Fixed-point
  list(APPEND unit_tests "test_logging_arkode_mristep.cpp\;2 1 1") # ImEx Newton
                                                                   # + Dense
  list(APPEND unit_tests "test_logging_arkode_mristep.cpp\;2 1 0") # ImEx Newton
                                                                   # + GMRES
  list(APPEND unit_tests "test_logging_arkode_mristep.cpp\;2 0") # ImEx
                                                                 # Fixed-point
  # MRIStep -- MRI-SR
  list(APPEND unit_tests "test_logging_arkode_mristep.cpp\;3") # Explicit
  list(APPEND unit_tests "test_logging_arkode_mristep.cpp\;4 1 1") # Implicit
                                                                   # Newton +
                                                                   # Dense
  list(APPEND unit_tests "test_logging_arkode_mristep.cpp\;4 1 0") # Implicit
                                                                   # Newton +
                                                                   # GMRES
  list(APPEND unit_tests "test_logging_arkode_mristep.cpp\;4 0") # Implicit
                                                                 # Fixed-point
  list(APPEND unit_tests "test_logging_arkode_mristep.cpp\;5 1 1") # ImEx Newton
                                                                   # + Dense
  list(APPEND unit_tests "test_logging_arkode_mristep.cpp\;5 1 0") # ImEx Newton
                                                                   # + GMRES
  list(APPEND unit_tests "test_logging_arkode_mristep.cpp\;5 0") # ImEx
                                                                 # Fixed-point
  # MRIStep -- MERK
  list(APPEND unit_tests "test_logging_arkode_mristep.cpp\;6") # Explicit
  # SplittingStep
  list(APPEND unit_tests "test_logging_arkode_splittingstep.cpp\;")
  # SPRKStep
  list(APPEND unit_tests "test_logging_arkode_sprkstep.cpp\;0")
  list(APPEND unit_tests "test_logging_arkode_sprkstep.cpp\;1") # Compensated
                                                                # sum
endif()

if(BUILD_CVODE)
  list(APPEND unit_tests "test_logging_cvode.cpp\;1 1") # Newton + Dense
  list(APPEND unit_tests "test_logging_cvode.cpp\;1 0") # Newton + GMRES
  list(APPEND unit_tests "test_logging_cvode.cpp\;0") # Fixed-point
endif()

if(BUILD_CVODES)
  list(APPEND unit_tests "test_logging_cvodes.cpp\;1 1") # Newton + Dense
  list(APPEND unit_tests "test_logging_cvodes.cpp\;1 0") # Newton + GMRES
  list(APPEND unit_tests "test_logging_cvodes.cpp\;0") # Fixed-point
endif()

if(BUILD_IDA)
  list(APPEND unit_tests "test_logging_ida.cpp\;1") # Newton + Dense
  list(APPEND unit_tests "test_logging_ida.cpp\;0") # Newton + GMRES
endif()

if(BUILD_IDAS)
  list(APPEND unit_tests "test_logging_idas.cpp\;1") # Newton + Dense
  list(APPEND unit_tests "test_logging_idas.cpp\;0") # Newton + GMRES
endif()

# if(BUILD_KINSOL) list(APPEND unit_tests "test_logging_kinsol.cpp\;") endif()

# Add the build and install targets for each test
foreach(test_tuple ${unit_tests})

  # parse the test tuple
  list(GET test_tuple 0 test_src)
  list(GET test_tuple 1 test_args)

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

  # Set the test name
  if(${SUNDIALS_LOGGING_LEVEL} GREATER 2)
    set(test_name "${test_target}_lvl${SUNDIALS_LOGGING_LEVEL}")
  else()
    set(test_name "${test_target}")
  endif()

  if("${test_args}" STREQUAL "")
    set(test_name ${test_name})
  else()
    string(REPLACE " " "_" test_name "${test_name}_${test_args}")
    string(REPLACE " " ";" test_args "${test_args}")
  endif()

  # add test to regression tests
  sundials_add_test(
    ${test_name} ${test_target}
    TEST_ARGS ${test_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    LABELS "logging")

endforeach()

message(STATUS "Added logging units tests")
