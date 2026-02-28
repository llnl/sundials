# ---------------------------------------------------------------
# Programmer(s):  Daniel R. Reynolds @ UMBC
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

# Add the build targets for each backend
if(BUILD_ARKODE
   AND BUILD_CVODE
   AND BUILD_IDA)
  foreach(backend ${KOKKOS_EXAMPLES_BACKENDS})

    # set benchmark target name
    set(benchmark_target "advection_reaction_3D_kokkos.${backend}")

    # benchmark source files
    sundials_add_executable(
      ${benchmark_target}
      advection_reaction_3D.cpp
      arkode_driver.cpp
      cvode_driver.cpp
      ida_driver.cpp
      rhs3D.hpp
      ParallelGrid.hpp
      check_retval.h)

    # which backend to use
    target_compile_definitions(${benchmark_target} PRIVATE USE_${backend})

    # directories to include
    target_include_directories(
      ${benchmark_target} PRIVATE ${PROJECT_SOURCE_DIR}/utilities
                                  ${MPI_CXX_INCLUDE_DIRS})

    # libraries to link against
    target_link_libraries(
      ${benchmark_target}
      PRIVATE sundials_arkode
              sundials_cvode
              sundials_ida
              sundials_nvecmpiplusx
              sundials_nveckokkos
              ${MPI_CXX_LIBRARIES}
              ${EXE_EXTRA_LINK_LIBS})

    install(
      TARGETS ${benchmark_target}
      DESTINATION "${BENCHMARKS_INSTALL_PATH}/advection_reaction_3D/kokkos")

    install(
      FILES README.md ../scripts/compare_error.py ../scripts/compute_error.py
            ../scripts/pickle_solution_output.py
      DESTINATION "${BENCHMARKS_INSTALL_PATH}/advection_reaction_3D/kokkos")

  endforeach()
endif()
