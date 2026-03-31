# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025, Lawrence Livermore National Security,
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

# test source files
set(target "ark_shampine")
sundials_add_executable(${target} ${target}.cpp)
set_target_properties(${target} PROPERTIES FOLDER "${target}")
target_include_directories(${target} PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/test)
target_link_libraries(${target} sundials_arkode sundials_nvecserial ${EXE_EXTRA_LINK_LIBS})

set(target "ark_vdp")
sundials_add_executable(${target} ${target}.cpp)
set_target_properties(${target} PROPERTIES FOLDER "${target}")
target_include_directories(${target} PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/test)
target_link_libraries(${target} sundials_arkode sundials_nvecserial ${EXE_EXTRA_LINK_LIBS})

set(target "ark_rlc_ode")
sundials_add_executable(${target} ${target}.cpp)
set_target_properties(${target} PROPERTIES FOLDER "${target}")
target_include_directories(${target}
  PRIVATE
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_BINARY_DIR}/src # need for sundials_lapack_def.h
  ${CMAKE_SOURCE_DIR}/test)
set(_lapack_targets LAPACK::LAPACK)
if(CMAKE_VERSION VERSION_LESS 3.20)
   list(APPEND _lapack_targets BLAS::BLAS)
endif()
target_link_libraries(${target} sundials_arkode sundials_nvecserial ${_lapack_targets} ${EXE_EXTRA_LINK_LIBS})

# Create symlink to plot.py in the build directory
file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/plot_shampine.py" "${CMAKE_CURRENT_BINARY_DIR}/plot_shampine.py" SYMBOLIC)
file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/plot_vdp.py" "${CMAKE_CURRENT_BINARY_DIR}/plot_vdp.py" SYMBOLIC)
file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/driver.sh" "${CMAKE_CURRENT_BINARY_DIR}/driver.sh" SYMBOLIC)

file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/plot_state.py" "${CMAKE_CURRENT_BINARY_DIR}/plot_state.py" SYMBOLIC)

# # add test to regression tests
# sundials_add_test(
#   ${test_name} ${test_target}
#   TEST_ARGS ${test_args}
#   ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
#   ANSWER_FILE ${test_name}.out
#   LABELS "logging")

message(STATUS "Added stiffness detection tests")
