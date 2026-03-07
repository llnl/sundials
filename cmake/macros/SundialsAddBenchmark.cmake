# ------------------------------------------------------------------------------
# Programmer(s): Yu Pan @ LLNL
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
# sundials_add_benchmark(<benchmark name> <executable> <base benchmark name>
#                        [NUM_CORES cores]
#                        [BENCHMARK_ARGS args]
#                        [IDENTIFIER ident]
#                        [TEST_RUNNER_ARGS args]
#                        [ENABLE_GPU]
# ~~~
#
# CMake function for adding benchmarks to `make benchmark`.
#
# The option NUM_CORES sets the number of cores (GPU count or CPU count) to run
# on/number of resource sets.
#
# The option BENCHMARK_ARGS sets the arguments to pass to the executable.
#
# The option IDENTIFIER sets the suffix to append to end of benchmark name
#
# The option TEST_RUNNER_ARGS sets the command line arguments to pass to the
# test executable
#
# The option ENABLE_GPU indicates this benchmark should be run with GPUs

function(sundials_add_benchmark NAME EXECUTABLE BASE_BENCHMARK_NAME)

  set(oneValueArgs NUM_CORES BENCHMARK_ARGS IDENTIFIER)
  set(multiValueArgs TEST_RUNNER_ARGS)
  set(options ENABLE_GPU)

  # parse inputs and create variables arg_<keyword>
  cmake_parse_arguments(arg "${options}" "${oneValueArgs}" "${multiValueArgs}"
                        ${ARGN})

  # set the target name
  if(arg_IDENTIFIER)
    set(TARGET_NAME ${NAME}_${arg_IDENTIFIER})
  else()
    if(arg_BENCHMARK_ARGS)
      string(REPLACE " " "_" TEST_SUFFIX "${arg_BENCHMARK_ARGS}")
      set(TARGET_NAME ${NAME}_${TEST_SUFFIX})
    else()
      set(TARGET_NAME ${NAME}_run)
    endif()
  endif()

  # make the the output directory if it doesn't exist
  set(_output_dir
      "${SUNDIALS_BENCHMARK_OUTPUT_DIR}/${BASE_BENCHMARK_NAME}/${TARGET_NAME}")
  if(NOT EXISTS ${_output_dir})
    file(MAKE_DIRECTORY ${_output_dir})
  endif()

  # command line arguments for the test runner script
  set(TEST_RUNNER_ARGS
      "--verbose" "--executablename=$<TARGET_FILE:${EXECUTABLE}>"
      "--outputdir=${_output_dir}" "--nodiff")

  if(SUNDIALS_ENABLE_CALIPER)
    set(_caliper_dir
        "${SUNDIALS_BENCHMARK_CALIPER_OUTPUT_DIR}/${BASE_BENCHMARK_NAME}/${TARGET_NAME}"
    )
    if(NOT EXISTS ${_caliper_dir})
      file(MAKE_DIRECTORY ${_caliper_dir})
    endif()
    list(APPEND TEST_RUNNER_ARGS "--profile")
    list(APPEND TEST_RUNNER_ARGS "--calidir=${_caliper_dir}")
  endif()

  # incorporate scheduler arguments into test_runner
  if(SUNDIALS_SCHEDULER_COMMAND STREQUAL "flux run")
    set(SCHEDULER_STRING " -n${arg_NUM_CORES}")
  elseif(SUNDIALS_SCHEDULER_COMMAND STREQUAL "jsrun")
    if(${arg_ENABLE_GPU})
      set(SCHEDULER_STRING
          " --smpiargs=\\\"-gpu\\\" -n${arg_NUM_CORES} -a1 -c1 -g1")
    else()
      set(SCHEDULER_STRING " -n${arg_NUM_CORES} -a1 -c1")
    endif()
  elseif(SUNDIALS_SCHEDULER_COMMAND STREQUAL "srun")
    set(SCHEDULER_STRING
        " -n${arg_NUM_CORES} --cpus-per-task=1 --ntasks-per-node=1")
  endif()
  string(REPLACE " " ";" SCHEDULER_ARGS "${SCHEDULER_STRING}")
  string(REPLACE " " ";" SCHEDULER_COMMAND_ARGS "${SUNDIALS_SCHEDULER_COMMAND}")

  string(STRIP "${RUN_COMMAND}" RUN_COMMAND)
  set(RUN_COMMAND ${SCHEDULER_COMMAND_ARGS} ${SCHEDULER_ARGS})
  list(APPEND TEST_RUNNER_ARGS "--runcommand=\"${RUN_COMMAND}\"")

  list(APPEND TEST_RUNNER_ARGS "--runargs=${arg_BENCHMARK_ARGS}"
       "--testname=${TARGET_NAME}")
  add_custom_target(
    ${TARGET_NAME}
    COMMENT "Running ${TARGET_NAME}"
    COMMAND ${PYTHON_EXECUTABLE} ${TESTRUNNER} ${TEST_RUNNER_ARGS})
  add_dependencies(benchmark ${TARGET_NAME})

endfunction()
