# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
#                modified by Daniel M. Margolis @ UMBC
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
# CMakeLists.txt file for the ARKode F2003 module serial examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;type" where the type is 'develop' for examples
# excluded from 'make test' in releases

set(FARKODE_examples "ark_analytic_f2003\;\;develop"
                     "ark_kpr_nestedmri_f2003\;\;develop")

# Regression tests
set(FARKODE_tests "test_ark_butcher_f2003\;develop")

if(SUNDIALS_INDEX_SIZE MATCHES "64")

  # Examples using SUNDIALS linear solvers
  list(
    APPEND
    FARKODE_examples
    "ark_bruss_f2003\;\;develop"
    "ark_diurnal_kry_bp_f2003\;\;develop"
    "ark_roberts_dns_f2003\;\;develop"
    "ark_kpr_mri_f2003\;\;develop"
    "ark_kpr_mri_f2003\;0 0.002\;develop"
    "ark_kpr_mri_f2003\;1 0.002\;develop"
    "ark_kpr_mri_f2003\;2 0.005\;develop"
    "ark_kpr_mri_f2003\;3 0.01\;develop"
    "ark_kpr_mri_f2003\;4 0.002\;develop"
    "ark_kpr_mri_f2003\;5 0.002\;develop"
    "ark_kpr_mri_f2003\;6 0.005\;develop"
    "ark_kpr_mri_f2003\;7 0.001\;develop"
    "ark_kpr_mri_f2003\;8 0.001\;develop"
    "ark_kpr_mri_f2003\;9 0.001\;develop")

  set(FARKODE_examples_KLU "ark_bruss1D_FEM_klu_f2003\;develop")

  set(FARKODE_examples_LAPACK "ark_roberts_dnsL_f2003\;\;develop")

endif()

# Specify libraries to link against
set(ARKODE_LIB sundials_arkode)
set(FARKODE_LIB sundials_farkode_mod)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${FARKODE_LIB} ${ARKODE_LIB} ${EXE_EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${FARKODE_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  if(NOT TARGET ${example})
    # Install fortran modules to a unique directory to avoid naming collisions
    set(CMAKE_Fortran_MODULE_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${example}.dir)

    # example source files
    sundials_add_executable(${example} ${example}.f90)

    # folder for IDEs
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS})
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
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE ${example_type})

  # install example source and out files
  if(SUNDIALS_EXAMPLES_ENABLE_INSTALL)
    install(FILES ${example}.f90 ${test_name}.out
            DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/arkode/F2003_serial)
  endif()
endforeach(example_tuple ${FARKODE_examples})

# Add the build and install targets for each KLU example (if needed)
if(BUILD_SUNLINSOL_KLU)

  # Sundials KLU linear solver module
  set(SUNLINSOLKLU_LIBS sundials_sunlinsolklu sundials_fsunlinsolklu_mod)

  # KLU libraries
  list(APPEND SUNLINSOLKLU_LIBS)

  foreach(example_tuple ${FARKODE_examples_KLU})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 example_type)

    # build fortran modules into a unique directory to avoid naming collisions
    set(CMAKE_Fortran_MODULE_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${example}.dir)

    # create the cmake executable target
    sundials_add_executable(${example} ${example}.f90)
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # add example to regression tests
    sundials_add_test(
      ${example} ${example}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${example}.out
      EXAMPLE_TYPE ${example_type})

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS} ${SUNLINSOLKLU_LIBS})

    # install example source and out files
    if(SUNDIALS_EXAMPLES_ENABLE_INSTALL)
      install(FILES ${example}.f90 ${example}.out
              DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/arkode/F2003_serial)
    endif()

  endforeach(example_tuple ${FARKODE_examples_KLU})

endif()

# Add the build and install targets for each LAPACK example (if needed)
if(BUILD_SUNLINSOL_LAPACKDENSE)

  # Sundials LAPACK linear solver modules
  set(SUNLINSOLLAPACK_LIBS sundials_sunlinsollapackdense
                           sundials_fsunlinsollapackdense_mod)

  foreach(example_tuple ${FARKODE_examples_LAPACK})

    # parse the example tuple
    list(GET example_tuple 0 example)
    list(GET example_tuple 1 example_type)

    # build fortran modules into a unique directory to avoid naming collisions
    set(CMAKE_Fortran_MODULE_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${example}.dir)

    # example source files
    sundials_add_executable(${example} ${example}.f90)

    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # add example to regression tests
    sundials_add_test(
      ${example} ${example}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${example}.out
      EXAMPLE_TYPE ${example_type})

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS} ${SUNLINSOLLAPACK_LIBS})

    # install example source and out files
    if(SUNDIALS_EXAMPLES_ENABLE_INSTALL)
      install(FILES ${example}.f90 ${example}.out
              DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/arkode/F2003_serial)
    endif()

  endforeach(example_tuple ${FARKODE_examples_LAPACK})

endif()

# Add the build and install targets for regression test
foreach(example_tuple ${FARKODE_tests})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_type)

  # Install fortran modules to a unique directory to avoid naming collisions
  set(CMAKE_Fortran_MODULE_DIRECTORY
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${example}.dir)

  # example source files
  sundials_add_executable(${example} ${example}.f90)

  set_target_properties(${example} PROPERTIES FOLDER "Examples")

  # add example to regression tests
  sundials_add_test(
    ${example} ${example}
    TEST_ARGS ""
    EXAMPLE_TYPE ${example_type}
    NODIFF)

  # libraries to link against
  target_link_libraries(${example} ${SUNDIALS_LIBS})

  # install example source and out files
  if(SUNDIALS_EXAMPLES_ENABLE_INSTALL)
    install(FILES ${example}.f90
            DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/arkode/F2003_serial)
  endif()

endforeach(example_tuple ${FARKODE_tests})

# create Makefile and CMakeLists.txt for examples
if(SUNDIALS_EXAMPLES_ENABLE_INSTALL)

  # Install the extra files
  foreach(extrafile ${ARKODE_extras})
    install(FILES ${extrafile}
            DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/arkode/F2003_serial)
  endforeach()

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "ARKODE")
  set(SOLVER_LIB "sundials_arkode")
  set(SOLVER_FLIB "sundials_farkode_mod")

  examples2string(FARKODE_examples EXAMPLES)
  examples2string(FARKODE_tests TESTFILES)
  set(EXAMPLES "${EXAMPLES} ${TESTFILES}")

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then be used
  # as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_F2003_ex.in
    ${PROJECT_BINARY_DIR}/examples/arkode/F2003_serial/CMakeLists.txt @ONLY)

  # install CMakelists.txt
  install(
    FILES ${PROJECT_BINARY_DIR}/examples/arkode/F2003_serial/CMakeLists.txt
    DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/arkode/F2003_serial)

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template for the
  # user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_F2003_ex.in
      ${PROJECT_BINARY_DIR}/examples/arkode/F2003_serial/Makefile_ex @ONLY)
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/arkode/F2003_serial/Makefile_ex
      DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/arkode/F2003_serial
      RENAME Makefile)
  endif()

  # add test_install target
  sundials_add_test_install(arkode F2003_serial)

endif()
