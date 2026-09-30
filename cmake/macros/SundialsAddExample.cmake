# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# Add an example executable and, optionally, its primary test and install files.
#
# sundials_add_example(<name> <sources>...
#                      [ADD_TEST] [NODIFF]
#                      [LINK_LIBRARIES libraries...]
#                      [INCLUDE_DIRECTORIES directories...]
#                      [COMPILE_DEFINITIONS definitions...]
#                      [PROPERTIES properties...]
#                      [TEST_NAME name] [TEST_ARGS arguments...]
#                      [MPI_NPROCS count] [ANSWER_FILE file]
#                      [ANSWER_DIR directory] [EXAMPLE_TYPE type]
#                      [INSTALL] [INSTALL_FILES files...])
#
# The positional sources and the SCALAR_TYPE, LINK_LIBRARIES,
# INCLUDE_DIRECTORIES, COMPILE_DEFINITIONS, and PROPERTIES arguments are passed
# to sundials_add_executable.
#
# ADD_TEST adds the primary regression test for the executable. TEST_NAME
# overrides the executable name used for the test. TEST_ARGS, EXTRA_ARGS,
# MPI_NPROCS, FLOAT_PRECISION, INTEGER_PRECISION, ANSWER_DIR, ANSWER_FILE,
# EXAMPLE_TYPE, LABELS, and NODIFF are passed to sundials_add_test.
#
# When SUNDIALS_ENABLE_EXAMPLES_INSTALL is enabled, INSTALL installs the
# positional source files and INSTALL_FILES installs the positional sources
# plus any additional source, header, input, or answer files. Files are placed
# in the current example directory relative to PROJECT_SOURCE_DIR/examples.
# ---------------------------------------------------------------

function(sundials_add_example name)
  set(options ADD_TEST INSTALL NODIFF)
  set(one_value_args
      SCALAR_TYPE
      TEST_NAME
      MPI_NPROCS
      FLOAT_PRECISION
      INTEGER_PRECISION
      ANSWER_DIR
      ANSWER_FILE
      EXAMPLE_TYPE)
  set(multi_value_args
      LINK_LIBRARIES
      INCLUDE_DIRECTORIES
      COMPILE_DEFINITIONS
      PROPERTIES
      TEST_ARGS
      EXTRA_ARGS
      LABELS
      INSTALL_FILES)
  cmake_parse_arguments(arg "${options}" "${one_value_args}"
                        "${multi_value_args}" ${ARGN})

  set(_executable_args ${arg_UNPARSED_ARGUMENTS})
  foreach(_keyword LINK_LIBRARIES INCLUDE_DIRECTORIES COMPILE_DEFINITIONS
                   PROPERTIES)
    if(DEFINED arg_${_keyword} AND NOT "${arg_${_keyword}}" STREQUAL "")
      list(APPEND _executable_args ${_keyword} ${arg_${_keyword}})
    endif()
  endforeach()
  if(arg_SCALAR_TYPE)
    list(APPEND _executable_args SCALAR_TYPE ${arg_SCALAR_TYPE})
  endif()
  sundials_add_executable(${name} ${_executable_args})

  if(arg_ADD_TEST)
    if(arg_TEST_NAME)
      set(_test_name ${arg_TEST_NAME})
    else()
      set(_test_name ${name})
    endif()
    set(_test_args)
    foreach(_keyword MPI_NPROCS FLOAT_PRECISION INTEGER_PRECISION ANSWER_DIR
                     ANSWER_FILE EXAMPLE_TYPE)
      if(DEFINED arg_${_keyword} AND NOT "${arg_${_keyword}}" STREQUAL "")
        list(APPEND _test_args ${_keyword} ${arg_${_keyword}})
      endif()
    endforeach()
    foreach(_keyword TEST_ARGS EXTRA_ARGS LABELS)
      if(DEFINED arg_${_keyword} AND NOT "${arg_${_keyword}}" STREQUAL "")
        list(APPEND _test_args ${_keyword} ${arg_${_keyword}})
      endif()
    endforeach()
    if(arg_NODIFF)
      list(APPEND _test_args NODIFF)
    endif()
    sundials_add_test(${_test_name} ${name} ${_test_args})
  endif()

  if(SUNDIALS_ENABLE_EXAMPLES_INSTALL AND (arg_INSTALL OR arg_INSTALL_FILES))
    file(RELATIVE_PATH _install_subdir "${PROJECT_SOURCE_DIR}/examples"
         "${CMAKE_CURRENT_SOURCE_DIR}")
    set(_install_files ${arg_UNPARSED_ARGUMENTS} ${arg_INSTALL_FILES})
    list(REMOVE_DUPLICATES _install_files)
    install(FILES ${_install_files}
            DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/${_install_subdir})
  endif()
endfunction()
