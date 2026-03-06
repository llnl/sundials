# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
# CMakeLists.txt file for the IDA F2003 module parallel examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;nodes\;tasks\;type" where the type is develop
# for examples excluded from 'make test' in releases
if(SUNDIALS_INDEX_SIZE MATCHES "64")
  set(FIDA_examples "ida_heat2D_kry_bbd_f2003\;\;1\;4\;develop\;2")
endif()

# Specify libraries to link against
set(IDA_LIBS sundials_ida sundials_fida_mod sundials_nvecparallel
             sundials_fnvecparallel_mod)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${IDA_LIBS} ${EXE_EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${FIDA_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 number_of_nodes)
  list(GET example_tuple 3 number_of_tasks)
  list(GET example_tuple 4 example_type)
  list(GET example_tuple 5 example_precision)

  if(NOT TARGET ${example})
    # Install fortran modules to a unique directory to avoid naming collisions
    set(CMAKE_Fortran_MODULE_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${example}.dir)

    # example source files
    sundials_add_executable(${example} ${example}.f90)

    # libraries to link against
    target_link_libraries(${example} MPI::MPI_Fortran ${SUNDIALS_LIBS})

    set_target_properties(${example} PROPERTIES FOLDER "Examples")
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
    MPI_NPROCS ${number_of_tasks}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE ${example_type}
    FLOAT_PRECISION ${example_precision})

  # find all .out files for this example
  file(GLOB example_out ${example}*.out)

  # install example source and out files
  if(SUNDIALS_EXAMPLES_ENABLE_INSTALL)
    install(FILES ${example}.f90 ${example_out}
            DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/ida/F2003_parallel)
  endif()

endforeach(example_tuple ${FIDA_examples})

# create Makefile and CMakeLists.txt for examples
if(SUNDIALS_EXAMPLES_ENABLE_INSTALL)

  # Install the extra files
  foreach(extrafile ${IDA_extras})
    install(FILES ${extrafile}
            DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/ida/F2003_parallel)
  endforeach()

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "IDA")

  # Makefile: convert semi-colon separated target list to space separated string
  list2string(IDA_LIBS EXAMPLE_LIBS)

  # CMakeLists: replace sundials_ prefix and convert to space separated string
  list(TRANSFORM IDA_LIBS REPLACE "sundials_" "SUNDIALS::"
                                  OUTPUT_VARIABLE EXAMPLES_CMAKE_TARGETS_tmp)
  list2string(EXAMPLES_CMAKE_TARGETS_tmp EXAMPLES_CMAKE_TARGETS)

  examples2string(FIDA_examples EXAMPLES)

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then be used
  # as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_parallel_F2003_ex.in
    ${PROJECT_BINARY_DIR}/examples/ida/F2003_parallel/CMakeLists.txt @ONLY)

  # install CMakelists.txt
  install(FILES ${PROJECT_BINARY_DIR}/examples/ida/F2003_parallel/CMakeLists.txt
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/ida/F2003_parallel)

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template for the
  # user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_parallel_F2003_ex.in
      ${PROJECT_BINARY_DIR}/examples/ida/F2003_parallel/Makefile_ex @ONLY)
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/ida/F2003_parallel/Makefile_ex
      DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/ida/F2003_parallel
      RENAME Makefile)
  endif()

  # add test_install target
  sundials_add_test_install(ida F2003_parallel)

endif()
