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
# Add an example executable, its primary test, and its install files.
#
# sundials_add_example(<name> SOURCES <sources>...
#                      [NO_TEST] [NO_INSTALL] [NO_DIFF]
#                      [LINK_LIBRARIES libraries...]
#                      [INCLUDE_DIRECTORIES directories...]
#                      [COMPILE_DEFINITIONS definitions...]
#                      [PROPERTIES properties...]
#                      [SOURCE_PROPERTIES properties...]
#                      [TEST_NAME name] [TEST_ARGS arguments...]
#                      [MPI_NPROCS count] [ANSWER_FILE file]
#                      [ANSWER_DIR directory] [EXAMPLE_TYPE type]
#                      [INSTALL_FILES files...])
#
# The SOURCES, SCALAR_TYPE, LINK_LIBRARIES,
# INCLUDE_DIRECTORIES, COMPILE_DEFINITIONS, PROPERTIES, and SOURCE_PROPERTIES
# arguments are passed to sundials_add_executable.
#
# A primary regression test is added unless NO_TEST is given. TEST_NAME defaults
# to the executable name, ANSWER_DIR to CMAKE_CURRENT_SOURCE_DIR, ANSWER_FILE to
# <test-name>.out, and EXAMPLE_TYPE to an empty value. TEST_ARGS, EXTRA_ARGS,
# MPI_NPROCS, FLOAT_PRECISION, INTEGER_PRECISION, LABELS, and NO_DIFF are passed
# to sundials_add_test.
#
# The executable's FOLDER property defaults to Examples. When
# SUNDIALS_ENABLE_EXAMPLES_INSTALL is enabled, the positional sources, primary
# answer file, and any additional source, header, input, or answer files listed
# in INSTALL_FILES are installed unless NO_INSTALL is given. Files are placed in
# the current example directory relative to PROJECT_SOURCE_DIR/examples.
# ---------------------------------------------------------------

function(sundials_add_example name)
  if(DEFINED SUNDIALS_EXAMPLE_VARIANT AND NOT name STREQUAL
                                          SUNDIALS_EXAMPLE_VARIANT)
    return()
  endif()

  set(options NO_INSTALL NO_TEST NO_DIFF)
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
      SOURCE_PROPERTIES
      TEST_ARGS
      EXTRA_ARGS
      LABELS
      INSTALL_FILES
      SOURCES)
  cmake_parse_arguments(arg "${options}" "${one_value_args}"
                        "${multi_value_args}" ${ARGN})

  if(arg_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "Unexpected arguments: ${arg_UNPARSED_ARGUMENTS}")
  endif()
  if(NOT arg_SOURCES)
    message(FATAL_ERROR "sundials_add_example(${name}) requires SOURCES")
  endif()

  set(_executable_args ${arg_SOURCES})
  foreach(_keyword LINK_LIBRARIES INCLUDE_DIRECTORIES COMPILE_DEFINITIONS
                   SOURCE_PROPERTIES)
    list(APPEND _executable_args ${_keyword} ${arg_${_keyword}})
  endforeach()
  list(APPEND _executable_args SCALAR_TYPE ${arg_SCALAR_TYPE})
  set(_properties ${arg_PROPERTIES})
  list(FIND _properties FOLDER _folder_index)
  if(_folder_index EQUAL -1)
    list(APPEND _properties FOLDER "Examples")
  endif()
  list(APPEND _executable_args PROPERTIES ${_properties})
  sundials_add_executable(${name} ${_executable_args})

  # A scalar-type mismatch intentionally does not create a target.
  if(NOT TARGET ${name})
    return()
  endif()

  if(arg_TEST_NAME)
    set(_test_name ${arg_TEST_NAME})
  else()
    set(_test_name ${name})
  endif()
  if(arg_ANSWER_FILE)
    set(_answer_file ${arg_ANSWER_FILE})
  else()
    set(_answer_file ${_test_name}.out)
  endif()

  if(NOT arg_NO_TEST)
    set(_test_args)
    foreach(_keyword MPI_NPROCS FLOAT_PRECISION INTEGER_PRECISION)
      list(APPEND _test_args ${_keyword} ${arg_${_keyword}})
    endforeach()
    if(arg_ANSWER_DIR)
      list(APPEND _test_args ANSWER_DIR ${arg_ANSWER_DIR})
    else()
      list(APPEND _test_args ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR})
    endif()
    list(APPEND _test_args ANSWER_FILE ${_answer_file})
    list(APPEND _test_args EXAMPLE_TYPE ${arg_EXAMPLE_TYPE})
    foreach(_keyword TEST_ARGS EXTRA_ARGS LABELS)
      list(APPEND _test_args ${_keyword} ${arg_${_keyword}})
    endforeach()
    if(arg_NO_DIFF)
      list(APPEND _test_args NO_DIFF)
    endif()
    sundials_add_test(${_test_name} ${name} ${_test_args})
  endif()

  if(SUNDIALS_ENABLE_EXAMPLES_INSTALL AND NOT arg_NO_INSTALL)
    file(RELATIVE_PATH _install_subdir "${PROJECT_SOURCE_DIR}/examples"
         "${CMAKE_CURRENT_SOURCE_DIR}")
    set(_install_files ${arg_SOURCES} ${arg_INSTALL_FILES})
    if(NOT arg_NO_TEST AND NOT arg_NO_DIFF)
      list(APPEND _install_files ${_answer_file})
    endif()
    list(REMOVE_DUPLICATES _install_files)
    install(FILES ${_install_files}
            DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/${_install_subdir})

    include(SundialsInstallExampleCMake)
    sundials_install_example_cmake(${name})
  endif()
endfunction()
