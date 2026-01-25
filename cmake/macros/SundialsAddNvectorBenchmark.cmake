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
# CMake macro for adding nvector benchmark executables.
# ---------------------------------------------------------------

function(sundials_add_nvector_benchmark NAME)

  set(options)
  set(singleValueArgs)
  set(multiValueArgs SOURCES SUNDIALS_TARGETS LINK_LIBRARIES INSTALL_SUBDIR)

  cmake_parse_arguments(arg "${options}" "${singleValueArgs}"
                        "${multiValueArgs}" ${ARGN})

  set(BENCHMARKS_DIR ${PROJECT_SOURCE_DIR}/benchmarks)

  add_executable(${NAME} ${BENCHMARKS_DIR}/nvector/test_nvector_performance.c
                         ${arg_SOURCES})

  set_target_properties(${NAME} PROPERTIES FOLDER "Benchmarks")

  target_include_directories(${NAME} PRIVATE ${BENCHMARKS_DIR}/nvector)

  target_link_libraries(${NAME} PRIVATE ${arg_SUNDIALS_TARGETS}
                                        ${arg_LINK_LIBRARIES} -lm)

  install(TARGETS ${NAME}
          DESTINATION "${BENCHMARKS_INSTALL_PATH}/${arg_INSTALL_SUBDIR}")

endfunction(sundials_add_nvector_benchmark)
