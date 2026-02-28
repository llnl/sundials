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

message(STATUS "Added OpenMPDEV NVECTOR benchmark")

set(NVECS_LIB sundials_nvecopenmpdev)

# Set-up linker flags and link libraries
set(SUNDIALS_LIBS ${NVECS_LIB} ${EXE_EXTRA_LINK_LIBS})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_C_FLAGS}")

sundials_add_executable(
  test_nvector_performance_openmpdev test_nvector_performance_openmpdev.c
  ../test_nvector_performance.c ../../../src/sundials/sundials_nvector.c)

# folder to organize targets in an IDE
set_target_properties(test_nvector_performance_openmp PROPERTIES FOLDER
                                                                 "Benchmarks")

target_include_directories(test_nvector_performance_openmpdev PRIVATE ..)

# libraries to link against
target_link_libraries(test_nvector_performance_openmpdev ${SUNDIALS_LIBS})

install(TARGETS test_nvector_performance_openmpdev
        DESTINATION "${BENCHMARKS_INSTALL_PATH}/nvector/openmpdev")
