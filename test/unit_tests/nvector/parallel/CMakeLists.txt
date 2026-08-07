# ---------------------------------------------------------------
# Programmer(s): Steven Smith @ LLNL
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
# CMakeLists.txt file for MPI nvector examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;nodes\;tasks\;type" where the type is
# develop for examples excluded from 'make test' in releases

# Examples using SUNDIALS MPI nvector
set(nvector_parallel_examples
    "test_nvector_mpi\;1000 0\;\;\;" # run sequentially
    "test_nvector_mpi\;1000 0\;1\;4\;" # run parallel on 4 procs
)

# Dependencies for nvector examples
set(nvector_examples_dependencies test_nvector test_mpinvector)

# If building F2003 tests
if(SUNDIALS_ENABLE_FORTRAN)
  set(nvector_parallel_fortran_examples "test_fnvector_parallel_mod\;\;1\;4\;")
endif()

# Add source directory to include directories
include_directories(. ..)

if(MPI_C_COMPILER)
  # use MPI wrapper as the compiler
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
else()
  # add MPI_INCLUDE_PATH to include directories
  include_directories(${MPI_INCLUDE_PATH})
endif()

if(SUNDIALS_ENABLE_FORTRAN)
  if(MPI_Fortran_COMPILER)
    # use MPI wrapper as the compiler
    set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
  else()
    # add MPI_INCLUDE_PATH to include directories
    include_directories(${MPI_INCLUDE_PATH})
  endif()
endif()

# Add the build and install targets for each example
foreach(example_tuple ${nvector_parallel_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 number_of_nodes)
  list(GET example_tuple 3 number_of_tasks)
  list(GET example_tuple 4 example_type)

  # check if this example has already been added, only need to add example
  # source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    sundials_add_executable(${example} ${example}.c)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(
      ${example} PRIVATE test_nvector_obj test_nvectormpi_obj
                         sundials_nvecparallel ${EXE_EXTRA_LINK_LIBS})

    if(NOT MPI_C_COMPILER)
      target_link_libraries(${example} PRIVATE ${MPI_LIBRARIES})
    endif()
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    if("${number_of_tasks}" STREQUAL "")
      set(test_name ${example})
    else()
      string(REGEX REPLACE " " "_" test_name ${example}_${number_of_tasks})
    endif()
  else()
    if("${number_of_tasks}" STREQUAL "")
      string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
    else()
      string(REGEX REPLACE " " "_" test_name
                           ${example}_${number_of_tasks}_${example_args})
    endif()
  endif()

  # add example to regression tests
  sundials_add_test(
    ${test_name} ${example}
    TEST_ARGS ${example_args}
    MPI_NPROCS ${number_of_tasks}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

endforeach(example_tuple ${nvector_parallel_examples})

# Add the build and install targets for each example
foreach(example_tuple ${nvector_parallel_fortran_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 number_of_nodes)
  list(GET example_tuple 3 number_of_tasks)
  list(GET example_tuple 4 example_type)

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
      ${SUNDIALS_SOURCE_DIR}/test/unit_tests/nvector/test_nvector.f90)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example} sundials_nvecparallel
                          sundials_fnvecparallel_mod ${EXE_EXTRA_LINK_LIBS})

    if(NOT MPI_Fortran_COMPILER)
      target_link_libraries(${example} ${MPI_LIBRARIES})
    endif()
  endif()

  # check if example args are provided and set the test name
  if("${example_args}" STREQUAL "")
    if("${number_of_tasks}" STREQUAL "")
      set(test_name ${example})
    else()
      string(REGEX REPLACE " " "_" test_name ${example}_${number_of_tasks})
    endif()
  else()
    if("${number_of_tasks}" STREQUAL "")
      string(REGEX REPLACE " " "_" test_name ${example}_${example_args})
    else()
      string(REGEX REPLACE " " "_" test_name
                           ${example}_${number_of_tasks}_${example_args})
    endif()
  endif()

  # add example to regression tests
  sundials_add_test(
    ${test_name} ${example}
    TEST_ARGS ${example_args}
    MPI_NPROCS ${number_of_tasks}
    EXAMPLE_TYPE ${example_type}
    NODIFF)

endforeach(example_tuple ${nvector_parallel_fortran_examples})
