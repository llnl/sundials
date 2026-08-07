# ---------------------------------------------------------------
# Programmer(s): Slaven Peles and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for Trilinos nvector examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;tasks\;type" where the type is 'develop'
# for examples excluded from 'make test' in releases

# Examples using SUNDIALS Trilinos nvector wrapper
if(MPI_CXX_FOUND)
  set(nvector_trilinos_examples
      "test_nvector_trilinos\;1000 0\;\;\;" # run sequentially
      "test_nvector_trilinos\;1000 0\;1\;4\;" # run parallel on 4 procs
  )
else()
  set(nvector_trilinos_examples
      "test_nvector_trilinos\;1000 0\;\;\;" # run sequentially
  )
endif()

# Dependencies for nvector examples
if(SUNDIALS_ENABLE_MPI)
  set(nvector_examples_dependencies test_nvector test_mpinvector)
else()
  set(nvector_examples_dependencies test_nvector)
endif()

# Add source directory to include directories
include_directories(. ..)

# Specify libraries to link against
set(NVECS_LIB sundials_nvectrilinos)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${NVECS_LIB} ${EXE_EXTRA_LINK_LIBS})

# Add the build and install targets for each example
foreach(example_tuple ${nvector_trilinos_examples})

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
    sundials_add_executable(${example} ${example}.cpp)

    # link vector test utilities
    target_link_libraries(${example} PRIVATE test_nvector_obj)
    if(SUNDIALS_ENABLE_MPI)
      target_link_libraries(${example} PRIVATE test_nvectormpi_obj)
    endif()

    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # Trilinos libraries must be linked after SUNDIALS libraries! Otherwise,
    # segfaults may occur; something to do with RCPs.
    target_link_libraries(${example} PRIVATE ${SUNDIALS_LIBS})
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

endforeach(example_tuple ${nvector_trilinos_examples})
