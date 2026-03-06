# ---------------------------------------------------------------
# Programmer(s): Ting Yan @ UMBC
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
# -----------------------------------------------------------------
# CMakeLists.txt file for IDAS OpenMP examples
# -----------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is 'develop' for
# examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(IDAS_examples "idasFoodWeb_bnd_omp\;4\;develop"
                  "idasFoodWeb_kry_omp\;4\;develop")

# Specify libraries to link against
set(IDAS_LIB sundials_idas)
set(NVECOMP_LIB sundials_nvecopenmp)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${IDAS_LIB} ${NVECOMP_LIB} ${EXE_EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${IDAS_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # check if this example has already been added, only need to add example
  # source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    sundials_add_executable(${example} ${example}.c)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example} ${SUNDIALS_LIBS})
  endif()

  # check if example args are provided and set the test name
  # IF("${example_args}" STREQUAL "") SET(test_name ${example}) ELSE()
  # STRING(REGEX REPLACE " " "_" test_name ${example}_${example_args}) ENDIF()

  # LEB NOTE: This goes away with TestRunner fix for OpenMP
  if(NOT ${example} MATCHES "idasFoodWeb_kry_omp")
    # add example to regression tests
    sundials_add_test(
      ${example} ${example}
      TEST_ARGS ${example_args}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ${example}.out
      EXAMPLE_TYPE ${example_type})
  endif()

  # install example source and .out files
  if(SUNDIALS_EXAMPLES_ENABLE_INSTALL)
    install(FILES ${example}.c ${example}.out
            DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/idas/C_openmp)
  endif()

endforeach(example_tuple ${IDAS_examples})

# create Makefile and CMakeLists.txt for examples
if(SUNDIALS_EXAMPLES_ENABLE_INSTALL)

  # Install the README file
  install(FILES README
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/idas/C_openmp)

  # Install the extra files
  foreach(extrafile ${IDAS_extras})
    install(FILES ${extrafile}
            DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/idas/C_openmp)
  endforeach()

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "IDAS")
  set(SOLVER_LIB "sundials_idas")

  examples2string(IDAS_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then be used
  # as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_openmp_C_ex.in
    ${PROJECT_BINARY_DIR}/examples/idas/C_openmp/CMakeLists.txt @ONLY)

  # install CMakelists.txt
  install(FILES ${PROJECT_BINARY_DIR}/examples/idas/C_openmp/CMakeLists.txt
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/idas/C_openmp)

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template for the
  # user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_openmp_C_ex.in
      ${PROJECT_BINARY_DIR}/examples/idas/C_openmp/Makefile_ex @ONLY)
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/idas/C_openmp/Makefile_ex
      DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/idas/C_openmp
      RENAME Makefile)
  endif()

  # add test_install target
  sundials_add_test_install(idas C_openmp)

endif()
