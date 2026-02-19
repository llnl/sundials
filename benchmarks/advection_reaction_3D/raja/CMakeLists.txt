# ---------------------------------------------------------------
# Programmer(s):  Cody J. Balos @ LLNL
#                 Daniel R. Reynolds @ UMBC
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

if(BUILD_ARKODE
   AND BUILD_CVODE
   AND BUILD_IDA)

  if((RAJA_BACKENDS MATCHES "TARGET_OPENMP") OR (RAJA_BACKENDS MATCHES "OPENMP"
                                                ))
    set(OTHER_LIBS OpenMP::OpenMP_CXX)
  endif()

  # Set up parameters to run benchmarks with
  set(BENCHMARK_VAR
      "--method ARK-IMEX --nls tl-newton --tf 0.01 --dont-save\;arkimex_tlnewton"
      "--method ARK-DIRK --nls newton --tf 0.01 --dont-save\;arkdirk_newton"
      "--method CV-BDF --nls newton --tf 0.01 --dont-save\;cvbdf_newton"
      "--method IDA --nls newton --tf 0.01 --dont-save\;ida_newton")

  # ----------------------------------------------------------------------------
  # MPI only
  # ----------------------------------------------------------------------------

  sundials_add_executable(
    advection_reaction_3D_raja
    advection_reaction_3D.cpp
    arkode_driver.cpp
    cvode_driver.cpp
    ida_driver.cpp
    rhs3D.hpp
    ParallelGrid.hpp
    check_retval.h
    backends.hpp)

  # ensure the linker language is reset to CXX
  set_target_properties(advection_reaction_3D_raja PROPERTIES LINKER_LANGUAGE
                                                              CXX)

  target_include_directories(
    advection_reaction_3D_raja PRIVATE ${PROJECT_SOURCE_DIR}/utilities
                                       ${MPI_CXX_INCLUDE_DIRS})

  target_link_libraries(
    advection_reaction_3D_raja
    PRIVATE sundials_arkode
            sundials_cvode
            sundials_ida
            sundials_nvecmpiplusx
            sundials_nvecserial
            RAJA
            ${MPI_CXX_LIBRARIES}
            ${OTHER_LIBS})

  install(TARGETS advection_reaction_3D_raja
          DESTINATION "${BENCHMARKS_INSTALL_PATH}/advection_reaction_3D/raja")

  install(FILES README.md ../scripts/compare_error.py
                ../scripts/compute_error.py ../scripts/pickle_solution_output.py
          DESTINATION "${BENCHMARKS_INSTALL_PATH}/advection_reaction_3D/raja")

  foreach(benchmark_tuple ${BENCHMARK_VAR})
    list(GET benchmark_tuple 0 benchmark_args)
    list(GET benchmark_tuple 1 identifier)

    sundials_add_benchmark(
      advection_reaction_3D_raja advection_reaction_3D_raja
      advection_reaction_3D
      NUM_CORES ${SUNDIALS_BENCHMARK_NUM_CPUS}
      BENCHMARK_ARGS ${benchmark_args}
      IDENTIFIER ${identifier})
  endforeach()

  # ----------------------------------------------------------------------------
  # MPI + CUDA
  # ----------------------------------------------------------------------------

  if(BUILD_NVECTOR_CUDA)

    set_source_files_properties(advection_reaction_3D.cpp PROPERTIES LANGUAGE
                                                                     CUDA)
    set_source_files_properties(arkode_driver.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(cvode_driver.cpp PROPERTIES LANGUAGE CUDA)
    set_source_files_properties(ida_driver.cpp PROPERTIES LANGUAGE CUDA)

    sundials_add_executable(
      advection_reaction_3D_raja_mpicuda
      advection_reaction_3D.cpp
      arkode_driver.cpp
      cvode_driver.cpp
      ida_driver.cpp
      rhs3D.hpp
      ParallelGrid.hpp
      check_retval.h
      backends.hpp)

    # ensure the linker language is reset to CXX
    set_target_properties(advection_reaction_3D_raja_mpicuda
                          PROPERTIES LINKER_LANGUAGE CXX)

    target_include_directories(
      advection_reaction_3D_raja_mpicuda PRIVATE ${PROJECT_SOURCE_DIR}/utilities
                                                 ${MPI_CXX_INCLUDE_DIRS})

    target_link_libraries(
      advection_reaction_3D_raja_mpicuda
      PRIVATE sundials_arkode
              sundials_cvode
              sundials_ida
              sundials_nvecmpiplusx
              sundials_nveccuda
              RAJA
              ${MPI_CXX_LIBRARIES}
              ${OTHER_LIBS})

    target_compile_definitions(advection_reaction_3D_raja_mpicuda
                               PRIVATE USE_CUDA_NVEC)

    install(TARGETS advection_reaction_3D_raja_mpicuda
            DESTINATION "${BENCHMARKS_INSTALL_PATH}/advection_reaction_3D/raja")

    foreach(benchmark_tuple ${BENCHMARK_VAR})
      list(GET benchmark_tuple 0 benchmark_args)
      list(GET benchmark_tuple 1 identifier)

      sundials_add_benchmark(
        advection_reaction_3D_raja_mpicuda advection_reaction_3D_raja_mpicuda
        advection_reaction_3D ENABLE_GPU
        NUM_CORES ${SUNDIALS_BENCHMARK_NUM_GPUS}
        BENCHMARK_ARGS ${benchmark_args}
        IDENTIFIER ${identifier})
    endforeach()
  endif()

  # ----------------------------------------------------------------------------
  # MPI + HIP
  # ----------------------------------------------------------------------------

  if(BUILD_NVECTOR_HIP)

    sundials_add_executable(
      advection_reaction_3D_raja_mpihip
      advection_reaction_3D.cpp
      advection_reaction_3D.hpp
      arkode_driver.cpp
      cvode_driver.cpp
      ida_driver.cpp
      rhs3D.hpp
      ParallelGrid.hpp
      check_retval.h
      backends.hpp)

    target_include_directories(
      advection_reaction_3D_raja_mpihip PRIVATE ${PROJECT_SOURCE_DIR}/utilities
                                                ${MPI_CXX_INCLUDE_DIRS})

    target_link_libraries(
      advection_reaction_3D_raja_mpihip
      PRIVATE sundials_arkode
              sundials_cvode
              sundials_ida
              sundials_nvecmpiplusx
              sundials_nvechip
              RAJA
              hip::device
              ${MPI_CXX_LIBRARIES}
              ${OTHER_LIBS})

    target_compile_definitions(advection_reaction_3D_raja_mpihip
                               PRIVATE USE_HIP_NVEC)

    install(TARGETS advection_reaction_3D_raja_mpihip
            DESTINATION "${BENCHMARKS_INSTALL_PATH}/advection_reaction_3D/raja")

    foreach(benchmark_tuple ${BENCHMARK_VAR})
      list(GET benchmark_tuple 0 benchmark_args)
      list(GET benchmark_tuple 1 identifier)

      sundials_add_benchmark(
        advection_reaction_3D_raja_mpihip advection_reaction_3D_raja_mpihip
        advection_reaction_3D ENABLE_GPU
        NUM_CORES ${SUNDIALS_BENCHMARK_NUM_GPUS}
        BENCHMARK_ARGS ${benchmark_args}
        IDENTIFIER ${identifier})
    endforeach()
  endif()
endif()
