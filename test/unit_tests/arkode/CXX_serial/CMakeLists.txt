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
# ARKODE C++ serial unit_tests
# ---------------------------------------------------------------

# List of test tuples of the form "name\;args\;type" where the type is 'develop'
# for tests excluded from 'make test' in releases
set(unit_tests
    "ark_test_accumerror_brusselator.cpp\;20 3 1\;exclude-single"
    "ark_test_accumerror_brusselator.cpp\;20 -4 0\;exclude-single"
    "ark_test_accumerror_brusselator.cpp\;20 5 0\;exclude-single"
    "ark_test_accumerror_kpr.cpp\;20 2 0\;exclude-single"
    "ark_test_accumerror_kpr.cpp\;20 3 1\;exclude-single"
    "ark_test_accumerror_kpr.cpp\;20 -4 1\;exclude-single"
    "ark_test_analytic_sys_mri.cpp\;0\;develop"
    "ark_test_analytic_sys_mri.cpp\;1\;develop"
    "ark_test_dahlquist_ark.cpp\;0 -1 0\;develop"
    "ark_test_dahlquist_ark.cpp\;0 0 0\;develop"
    "ark_test_dahlquist_ark.cpp\;0 0 1\;develop"
    "ark_test_dahlquist_ark.cpp\;0 1 0\;develop"
    "ark_test_dahlquist_ark.cpp\;0 1 1\;develop"
    "ark_test_dahlquist_ark.cpp\;1 -1 0\;develop"
    "ark_test_dahlquist_ark.cpp\;1 0 0\;develop"
    "ark_test_dahlquist_ark.cpp\;1 0 1\;develop"
    "ark_test_dahlquist_ark.cpp\;1 1 0\;develop"
    "ark_test_dahlquist_ark.cpp\;1 1 1\;develop"
    "ark_test_dahlquist_ark.cpp\;2 -1 0\;develop"
    "ark_test_dahlquist_ark.cpp\;2 0 0\;develop"
    "ark_test_dahlquist_ark.cpp\;2 0 1\;develop"
    "ark_test_dahlquist_ark.cpp\;2 1 0\;develop"
    "ark_test_dahlquist_ark.cpp\;2 1 1\;develop"
    "ark_test_dahlquist_erk.cpp\;-1\;develop"
    "ark_test_dahlquist_erk.cpp\;0\;develop"
    "ark_test_dahlquist_erk.cpp\;1\;develop"
    "ark_test_dahlquist_mri.cpp\;-1\;develop"
    "ark_test_dahlquist_mri.cpp\;0\;develop"
    "ark_test_dahlquist_mri.cpp\;1\;develop"
    "ark_test_butcher.cpp\;\;develop"
    "ark_test_getjac.cpp\;\;develop"
    "ark_test_getjac_mri.cpp\;\;develop"
    "ark_test_kpr_mriadapt.cpp\;--hs 0.002 --rtol 0.000004 --scontrol 0\;exclude-single"
    "ark_test_brusselator_mriadapt.cpp\;--rtol 0.000004 --scontrol 0\;exclude-single"
    "ark_test_slowerror_brusselator.cpp\;\;exclude-single"
    "ark_test_slowerror_kpr.cpp\;\;exclude-single"
    "ark_test_slowerror_polynomial.cpp\;\;exclude-single"
    "ark_test_splittingstep.cpp\;\;"
    "ark_test_adjoint_erk.cpp\;--check-freq 1\;"
    "ark_test_adjoint_erk.cpp\;--check-freq 2\;"
    "ark_test_adjoint_erk.cpp\;--check-freq 5\;"
    "ark_test_adjoint_erk.cpp\;--check-freq 1 --dont-keep\;"
    "ark_test_adjoint_erk.cpp\;--check-freq 2 --dont-keep\;"
    "ark_test_adjoint_erk.cpp\;--check-freq 5 --dont-keep\;"
    "ark_test_adjoint_ark.cpp\;--check-freq 1\;"
    "ark_test_adjoint_ark.cpp\;--check-freq 2\;"
    "ark_test_adjoint_ark.cpp\;--check-freq 5\;"
    "ark_test_adjoint_ark.cpp\;--check-freq 1 --dont-keep\;"
    "ark_test_adjoint_ark.cpp\;--check-freq 2 --dont-keep\;"
    "ark_test_adjoint_ark.cpp\;--check-freq 5 --dont-keep\;"
    "ark_test_prealloc_arkstep.cpp\;1\;"
    "ark_test_prealloc_erkstep.cpp\;1\;"
    "ark_test_prealloc_forcingstep.cpp\;1\;"
    "ark_test_prealloc_lsrkstep.cpp\;1\;"
    "ark_test_prealloc_mristep.cpp\;1\;"
    "ark_test_prealloc_splittingstep.cpp\;1\;"
    "ark_test_prealloc_sprkstep.cpp\;1\;")

# Add the build and install targets for each test
foreach(test_tuple ${unit_tests})

  # Parse the test tuple
  list(GET test_tuple 0 test)
  list(GET test_tuple 1 test_args)
  list(GET test_tuple 2 test_type)

  # Extract the file name without extension
  get_filename_component(test_target ${test} NAME_WE)

  # Check if this test has already been added, only need to add test source
  # files once for testing with different inputs
  if(NOT TARGET ${test_target})

    # Test source files
    sundials_add_executable(${test_target} ${test})

    # Folder to organize targets in an IDE
    set_target_properties(${test_target} PROPERTIES FOLDER "unit_tests")

    # Include location of public and private header files
    target_include_directories(
      ${test_target}
      PRIVATE $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
              ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/src
              ${CMAKE_SOURCE_DIR}/test/unit_tests)

    # We explicitly choose which object libraries to link to and link in the
    # arkode objects so that we have access to private functions w/o changing
    # their visibility in the installed libraries.
    target_link_libraries(
      ${test_target}
      PRIVATE $<TARGET_OBJECTS:sundials_arkode_obj>
              sundials_sunmemsys_obj
              sundials_nvecserial_obj
              sundials_nvecmanyvector_obj
              sundials_sunlinsolband_obj
              sundials_sunlinsoldense_obj
              sundials_sunlinsolspgmr_obj
              sundials_sunnonlinsolnewton_obj
              sundials_sunnonlinsolfixedpoint_obj
              sundials_sunadaptcontrollerimexgus_obj
              sundials_sunadaptcontrollersoderlind_obj
              sundials_sunadaptcontrollermrihtol_obj
              sundials_adjointcheckpointscheme_fixed_obj
              ${EXE_EXTRA_LINK_LIBS})

    # Tell CMake that we depend on the ARKODE library since it does not pick
    # that up from $<TARGET_OBJECTS:sundials_arkode_obj>.
    add_dependencies(${test_target} sundials_arkode_obj)

  endif()

  # Check if test args are provided and set the test name
  if("${test_args}" STREQUAL "")
    set(test_name ${test_target})
  else()
    string(REPLACE " " "_" test_name "${test_target}_${test_args}")
  endif()

  # add test to regression tests
  sundials_add_test(
    ${test_name} ${test_target}
    TEST_ARGS ${test_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE ${test_type})

endforeach()

message(STATUS "Added ARKODE CXX serial unit tests")
