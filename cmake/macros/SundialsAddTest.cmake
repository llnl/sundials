# ------------------------------------------------------------------------------
# Programmer(s): Steven Smith and David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------

# ~~~
# sundials_add_test(<test name> <executable>
#                   [NODIFF]
#                   [MPI_NPROCS num_processes]
#                   [FLOAT_PRECISION num_digits]
#                   [INTEGER_PRECISION percent_difference]
#                   [ANSWER_DIR path]
#                   [ANSWER_FIEL file]
#                   [EXAMPLE_TYPE type]
#                   [TEST_ARGS arg1 arg2 ...]
#                   [LABELS label1 label2 ...])
# ~~~
#
# CMake macro to add a SUNDIALS regression test. Keyword input arguments can be
# added after <executable> to set regression test options.
#
# The option NODIFF disables comparison of the test output against the answer
# file
#
# The option MPI_NPROCS sets the number of mpi tasks to use in parallel tests
#
# The option FLOAT_PRECISION set the precision (number of digits) for floating
# point failure comparisons. To use the default value, either don't provide the
# keyword, or provide the value "default".
#
# The option INTEGER_PRECISION sets the integer percentage difference for
# failure comparison.
#
# The option ANSWER_DIR sets the path to the directory containing the test
# answer file
#
# The option ANSWER_FILE set the name of test answer file
#
# The option EXAMPLE_TYPE set the example type i.e., release or develop examples
#
# The option TEST_ARGS are command line arguments to pass to the executable
#
# The options LABELS are labels added to the test properties to easily run (or
# exclude) groups of test with ctest -L <label> (or ctest -LE <label>)
#
# When SUNDIALS_TEST_ENABLE_DEV_TESTS is OFF (default) the executable is run and
# success or failure is determined by the executable return value (zero or
# non-zero respectively).
#
# When SUNDIALS_TEST_ENABLE_DEV_TESTS is ON the executable is run and its output
# is compared with the corresponding .out file. If the output differs
# significantly then the test fails. The default level of significance is 4
# decimal places for floating point values and 10% for integer values.
#
# The level of precision can be adjusted for an individual test with the
# FLOAT_PRECISION AND INTEGER_PRECISION keyword inputs to the macro or globally
# for all tests with the cache variables SUNDIALS_TEST_FLOAT_PRECISION and
# SUNDIALS_TEST_INTEGER_PRECISION.
#
# -D SUNDIALS_TEST_FLOAT_PRECISION=<number of digits>
#
# -D SUNDIALS_TEST_INTEGER_PRECISION=<% difference>
#
# By default testing output is written to builddir/Testing/output and the .out
# answer file directory is set using the ANSWER_DIR keyword input to
# sourcedir/examples/package/testdir. These can be changed by setting the cache
# variables SUNDIALS_TEST_OUTPUT_DIR and SUNDIALS_TEST_ANSWER_DIR.
#
# -D SUNDIALS_TEST_OUTPUT_DIR=<path to output directory>
#
# -D SUNDIALS_TEST_ANSWER_DIR=<path to answer directory>

function(SUNDIALS_ADD_TEST NAME EXECUTABLE)

  set(options "NODIFF")
  set(oneValueArgs "MPI_NPROCS" "FLOAT_PRECISION" "INTEGER_PRECISION"
                   "ANSWER_DIR" "ANSWER_FILE" "EXAMPLE_TYPE")
  set(multiValueArgs "LABELS" "TEST_ARGS" "EXTRA_ARGS")

  # parse inputs and create variables arg_<keyword>
  cmake_parse_arguments(arg "${options}" "${oneValueArgs}" "${multiValueArgs}"
                        ${ARGN})

  # ---------------------------------
  # check if the test should be added
  # ---------------------------------

  set(_add_test TRUE)

  # exclude development tests (non-empty example type)
  # TODO(DJG): When examples and development tests are separated this check can
  # be removed
  if(NOT SUNDIALS_TEST_ENABLE_DEV_TESTS AND arg_EXAMPLE_TYPE)
    set(_add_test FALSE)
  endif()

  # always excluded
  if("${arg_EXAMPLE_TYPE}" STREQUAL "exclude")
    set(_add_test FALSE)
  endif()

  # precision-specific exclusions
  string(TOLOWER "exclude-${SUNDIALS_PRECISION}" _exclude_precision)
  if("${arg_EXAMPLE_TYPE}" STREQUAL _exclude_precision)
    set(_add_test FALSE)
  endif()

  # --------
  # add test
  # --------

  if(_add_test)

    # ---------------------------------------
    # commands or flags before the executable
    # ---------------------------------------

    set(_pre_exe "")
    if(arg_MPI_NPROCS)
      if(SUNDIALS_TEST_MPIRUN_COMMAND)
        set(_pre_exe "${SUNDIALS_TEST_MPIRUN_COMMAND}")
      else()
        set(_pre_exe "${MPIEXEC_EXECUTABLE}")
      endif()
      if(MPIEXEC_PREFLAGS)
        set(_pre_exe "${_pre_exe} ${MPIEXEC_PREFLAGS}")
      endif()
      set(_pre_exe "${_pre_exe} ${MPIEXEC_NUMPROC_FLAG} ${arg_MPI_NPROCS}")
    endif()
    # Remove leading and trailing white space as it can cause erroneous test
    # failures with some MPI implementations
    string(STRIP "${_pre_exe}" _pre_exe)

    # ------------------------------------
    # flags or inputs after the executable
    # ------------------------------------

    # When checking if test command line args have been provided we compare
    # against an empty string (here and elsewhere below when checking _post_exe)
    # to avoid not adding the command line arg when it is also a false constant
    # in CMake (0, FALSE, OFF, etc.) e.g., "test_foo 0"
    set(_post_exe "")
    if(NOT "${arg_TEST_ARGS}" STREQUAL "")
      set(_post_exe "${arg_TEST_ARGS}")
    endif()
    if(NOT "${arg_EXTRA_ARGS}" STREQUAL "")
      if(NOT "${_post_exe}" STREQUAL "")
        set(_post_exe "${_post_exe} ${arg_EXTRA_ARGS}")
      else()
        set(_post_exe "${arg_EXTRA_ARGS}")
      endif()
    endif()
    if(arg_MPI_NPROCS AND MPIEXEC_POSTFLAGS)
      if(NOT "${_post_exe}" STREQUAL "")
        set(_post_exe "${MPIEXEC_POSTFLAGS} ${_post_exe}")
      else()
        set(_post_exe "${MPIEXEC_POSTFLAGS}")
      endif()
    endif()
    string(STRIP "${_post_exe}" _post_exe)

    # -------------------
    # create test command
    # -------------------

    # When using the test runner _pre_exe and _post_exe must be a string so we
    # replace semicolons with spaces below. Otherwise, _pre_exe and _post_exe
    # must be a semicolon separated list so we replace spaces with semicolons
    # below.

    if(SUNDIALS_TEST_USE_RUNNER)

      # command line arguments for the test runner script
      set(TEST_ARGS
          "--verbose" "--testname=${NAME}"
          "--executablename=$<TARGET_FILE:${EXECUTABLE}>"
          "--outputdir=${SUNDIALS_TEST_OUTPUT_DIR}")

      # set the test runcommand
      if(_pre_exe)
        string(REPLACE ";" " " _pre_exe "${_pre_exe}")
        list(APPEND TEST_ARGS "--runcommand=\"${_pre_exe}\"")
      endif()

      # set the test input args
      if(NOT "${_post_exe}" STREQUAL "")
        string(REPLACE ";" " " _post_exe "${_post_exe}")
        list(APPEND TEST_ARGS "--runargs=\"${_post_exe}\"")
      endif()

      if(SUNDIALS_TEST_ENABLE_PROFILING AND SUNDIALS_ENABLE_CALIPER)
        list(APPEND TEST_ARGS "--profile")
        list(APPEND TEST_ARGS "--calidir=${SUNDIALS_TEST_CALIPER_OUTPUT_DIR}")
      endif()

      # set comparison precisions or do not diff the output and answer files
      if(SUNDIALS_TEST_ENABLE_DIFF_OUTPUT AND NOT arg_NODIFF)

        # set answer directory
        if(SUNDIALS_TEST_ANSWER_DIR)
          list(APPEND TEST_ARGS "--answerdir=${SUNDIALS_TEST_ANSWER_DIR}")
        elseif(arg_ANSWER_DIR)
          list(APPEND TEST_ARGS "--answerdir=${arg_ANSWER_DIR}")
        endif()

        # set the test answer file name
        if(arg_ANSWER_FILE)
          list(APPEND TEST_ARGS "--answerfile=${arg_ANSWER_FILE}")
        endif()

        # set floating point precision
        if(arg_FLOAT_PRECISION AND (NOT arg_FLOAT_PRECISION MATCHES
                                    "DEFAULT|default"))
          list(APPEND TEST_ARGS "--floatprecision=${arg_FLOAT_PRECISION}")
        else()
          list(APPEND TEST_ARGS
               "--floatprecision=${SUNDIALS_TEST_FLOAT_PRECISION}")
        endif()

        # set integer precision
        if(arg_INTEGER_PRECISION AND (NOT arg_INTEGER_PRECISION MATCHES
                                      "DEFAULT|default"))
          list(APPEND TEST_ARGS "--integerpercentage=${arg_INTEGER_PRECISION}")
        else()
          list(APPEND TEST_ARGS
               "--integerpercentage=${SUNDIALS_TEST_INTEGER_PRECISION}")
        endif()

      else()

        list(APPEND TEST_ARGS "--nodiff")

      endif()

      add_test(NAME ${NAME} COMMAND ${Python3_EXECUTABLE} ${TESTRUNNER}
                                    ${TEST_ARGS})

    else()

      # set the test runcommand
      if(_pre_exe)
        string(REPLACE " " ";" _pre_exe "${_pre_exe}")
      endif()

      # set the test input args
      if(NOT "${_post_exe}" STREQUAL "")
        string(REPLACE " " ";" _post_exe "${_post_exe}")
      endif()

      add_test(NAME ${NAME} COMMAND ${_pre_exe} $<TARGET_FILE:${EXECUTABLE}>
                                    ${_post_exe})

    endif()

    # set any labels (must quote arg_LABELS)
    if(arg_LABELS)
      set_tests_properties(${NAME} PROPERTIES LABELS "${arg_LABELS}")
    endif()

  endif()

  unset(_add_test)
  unset(_pre_exe)
  unset(_post_exe)

endfunction()
