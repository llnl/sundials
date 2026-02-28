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

message(STATUS "Added Kokkos NVECTOR benchmark")

foreach(backend ${KOKKOS_EXAMPLES_BACKENDS})

  sundials_add_nvector_benchmark(
    test_nvector_performance_kokkos.${backend}
    SOURCES test_nvector_performance_kokkos.cpp
    SUNDIALS_TARGETS sundials_core sundials_nveckokkos
    INSTALL_SUBDIR nvector/kokkos)

  target_compile_definitions(test_nvector_performance_kokkos.${backend}
                             PRIVATE USE_${backend})

  install(TARGETS test_nvector_performance_kokkos.${backend}
          DESTINATION "${BENCHMARKS_INSTALL_PATH}/")

endforeach()
