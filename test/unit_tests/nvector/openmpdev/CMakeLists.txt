# ---------------------------------------------------------------
# Programmer(s): David J. Gardner and Shelby Lockhart @ LLNL
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
# CMakeLists.txt file for OpenMP DEV nvector examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is 'develop' for
# examples excluded from 'make test' in releases

# Examples using SUNDIALS OpenMP DEV nvector
set(nvector_openmpdev_examples "test_nvector_openmpdev\;1000 0\;"
                               "test_nvector_openmpdev\;10000 0\;")

# Dependencies for nvector examples
set(nvector_examples_dependencies test_nvector)

# Add source directory to include directories
include_directories(. ..)

# Specify libraries to link against
set(NVECS_LIB sundials_nvecopenmpdev)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${NVECS_LIB} ${EXE_EXTRA_LINK_LIBS})

set(CMAKE_EXE_LINKER_FLAGS ${OpenMP_C_FLAGS})

# Add the build and install targets for each example
foreach(example_tuple ${nvector_openmpdev_examples})

  # parse the example tuple
  list(GET example_tuple 0 example)
  list(GET example_tuple 1 example_args)
  list(GET example_tuple 2 example_type)

  # check if this example has already been added, only need to add example
  # source files once for testing with different inputs
  if(NOT TARGET ${example})
    # example source files
    sundials_add_executable(${example} ${example}.c)

    # link vector test utilities
    target_link_libraries(${example} PRIVATE test_nvector_obj)

    # folder to organize targets in an IDE
    set_target_properties(${example} PROPERTIES FOLDER "Examples")

    # libraries to link against
    target_link_libraries(${example} PRIVATE ${SUNDIALS_LIBS})
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

endforeach()
