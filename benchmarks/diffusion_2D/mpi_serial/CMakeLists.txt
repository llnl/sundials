# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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

# list of tests
set(tests)

if(BUILD_ARKODE)
  list(APPEND tests "arkode\;BENCHMARK_ODE")
endif()

if(BUILD_CVODE)
  list(APPEND tests "cvode\;BENCHMARK_ODE")
endif()

if(BUILD_IDA)
  list(APPEND tests "ida\;BENCHMARK_DAE")
endif()

# create executables
foreach(test_tuple ${tests})

  # parse the test tuple
  list(GET test_tuple 0 package)
  list(GET test_tuple 1 problem_type)

  set(sources ${benchmark_prefix}/main_${package}.cpp ${shared_sources}
              buffers.cpp diffusion.cpp solution.cpp utils.cpp)

  # set the target name
  set(target ${package}_diffusion_2D_mpi)

  # create executable
  sundials_add_executable(${target} ${sources})

  add_dependencies(benchmark ${target})

  set_target_properties(${target} PROPERTIES FOLDER "Benchmarks")

  target_compile_definitions(${target} PRIVATE ${problem_type})

  target_include_directories(${target} PRIVATE ${benchmark_prefix})

  target_link_libraries(${target} PRIVATE sundials_${package}
                                          sundials_nvecparallel MPI::MPI_CXX)

  if(BUILD_SUNLINSOL_SUPERLUDIST)
    target_compile_definitions(${target} PRIVATE USE_SUPERLU_DIST)
    target_link_libraries(${target} PRIVATE sundials_sunlinsolsuperludist
                                            sundials_sunmatrixslunrloc)
  endif()

  install(TARGETS ${target}
          DESTINATION "${BENCHMARKS_INSTALL_PATH}/diffusion_2D")

  sundials_add_benchmark(${target} ${target} diffusion_2D
                         NUM_CORES ${SUNDIALS_BENCHMARK_NUM_CPUS})

endforeach()
