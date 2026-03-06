# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ UMBC
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
# CMakeLists.txt file for sunlinsol SPFGMR examples
# ---------------------------------------------------------------

# Set tolerance for linear solver test based on Sundials precision
if(SUNDIALS_PRECISION MATCHES "SINGLE")
  set(TOL "1e-5")
elseif(SUNDIALS_PRECISION MATCHES "DOUBLE")
  set(TOL "1e-13")
else()
  set(TOL "1e-14")
endif()

# Example lists are tuples "name\;args\;type" where the type is 'develop' for
# examples excluded from 'make test' in releases

# Examples using SUNDIALS SPFGMR linear solver
set(sunlinsol_spfgmr_examples
    "test_sunlinsol_spfgmr_serial\;100 1 100 ${TOL} 0\;"
    "test_sunlinsol_spfgmr_serial\;100 2 100 ${TOL} 0\;")

# Dependencies for sunlinsol examples
set(sunlinsol_spfgmr_dependencies test_sunlinsol)

if(BUILD_FORTRAN_MODULE_INTERFACE)
  set(sunlinsol_spfgmr_fortran_examples "test_fsunlinsol_spfgmr_mod_serial\;\;")
endif()

# Add source directory to include directories
include_directories(. ../..)

# Add the build and install targets for each example
foreach(example_tuple ${sunlinsol_spfgmr_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # check if this example has already been added, only need to add example
  # source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    sundials_add_executable(${example} ${example}.c ../../test_sunlinsol.c)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example} sundials_nvecserial
                          sundials_sunlinsolspfgmr ${EXE_EXTRA_LINK_LIBS})
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    set(test_name ${example})
  else()
    string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
  endif()

  # add example to regression tests
  sundials_add_test(
    ${test_name} ${example}
    TEST_ARGS ${example_args}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

endforeach(example_tuple ${sunlinsol_spfgmr_examples})

# Add the build and install targets for each example
foreach(example_tuple ${sunlinsol_spfgmr_fortran_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # check if this example has already been added, only need to add example
  # source files once for testing with different inputs
  if(NOT TARGET ${example})
    # build fortran modules into a unique directory to avoid naming collisions
    set(CMAKE_Fortran_MODULE_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${example}.dir)

    # example source files
    sundials_add_executable(
      ${example} ${example}.f90
      ${SUNDIALS_SOURCE_DIR}/test/unit_tests/utilities/test_utilities.f90
      ${SUNDIALS_SOURCE_DIR}/test/unit_tests/sunlinsol/test_sunlinsol.f90)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(
      ${example} sundials_nvecserial sundials_fnvecserial_mod
      sundials_sunlinsolspfgmr sundials_fsunlinsolspfgmr_mod
      ${EXE_EXTRA_LINK_LIBS})

    # check if example args are provided and set the test name
    if("${example_args}" STREQUAL "")
      set(test_name ${example})
    else()
      string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
    endif()

    # add example to regression tests
    sundials_add_test(
      ${test_name} ${example}
      TEST_ARGS ${example_args}
      EXAMPLE_TYPE ${example_type}
      NODIFF)
  endif()

endforeach(example_tuple ${sunlinsol_spfgmr_fortran_examples})
