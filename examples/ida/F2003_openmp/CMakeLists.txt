# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ UMBC
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
# CMakeLists.txt file for the F2003 IDA OpenMP examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is 'develop' for
# examples excluded from 'make test' in releases

if(SUNDIALS_INDEX_SIZE MATCHES "64")
  # Examples using SUNDIALS linear solvers
  set(FIDA_examples_OMP "idaHeat2D_kry_omp_f2003.f90\;4\;exclude"
                        "idaHeat2D_kry_omp_f2003.f90\;8\;exclude")
endif()

# Specify libraries to link against
set(IDA_LIB sundials_fida_mod)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${IDA_LIB} ${EXE_EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${FIDA_examples_OMP})
  # parse the example tuple
  list(GET example_tuple 0 example_file)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # extract the file name without extension
  get_filename_component(example ${example_file} NAME_WE)

  if(NOT TARGET ${example})
    sundials_add_executable(${example} ${example_file})

    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(
      ${example} sundials_ida sundials_fida_mod sundials_nvecopenmp
      sundials_fnvecopenmp_mod ${SUNDIALS_LIBS})

    # Install fortran modules to a unique directory to avoid naming collisions
    set(CMAKE_Fortran_MODULE_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${example}.dir)

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
  if(EXAMPLES_INSTALL)
    install(FILES ${example}.f90 ${test_name}.out
            DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/F2003_openmp)
  endif()
endforeach(example_tuple ${FIDA_examples_OMP})

# create Makefile and CMakeLists.txt for examples
if(EXAMPLES_INSTALL)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "IDA")
  set(SOLVER_LIB "sundials_ida")
  set(SOLVER_FLIB "sundials_fida_mod")

  examples2string(FIDA_examples_OMP EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then be used
  # as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_openmp_F2003_ex.in
    ${PROJECT_BINARY_DIR}/examples/ida/F2003_openmp/CMakeLists.txt @ONLY)

  # install CMakelists.txt
  install(FILES ${PROJECT_BINARY_DIR}/examples/ida/F2003_openmp/CMakeLists.txt
          DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/F2003_openmp)

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template for the
  # user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_openmp_F2003_ex.in
      ${PROJECT_BINARY_DIR}/examples/ida/F2003_openmp/Makefile_ex @ONLY)
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/ida/F2003_openmp/Makefile_ex
      DESTINATION ${EXAMPLES_INSTALL_PATH}/ida/F2003_openmp
      RENAME Makefile)
  endif()

  # add test_install target
  sundials_add_test_install(ida F2003_openmp)

endif()
