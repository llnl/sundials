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

if(SUNDIALS_ENABLE_ARKODE
   OR SUNDIALS_ENABLE_CVODE
   OR SUNDIALS_ENABLE_IDA)

  # Shared sources
  set(shared_sources diffusion_2D.hpp diffusion_2D.cpp
                     preconditioner_jacobi.cpp)

  # Benchmark prefix
  set(benchmark_prefix ${SUNDIALS_SOURCE_DIR}/benchmarks/diffusion_2D/)

  add_prefix(${benchmark_prefix} shared_sources)

  # MPI + Serial
  add_subdirectory(mpi_serial)

  # MPI + CUDA/HIP
  if(SUNDIALS_ENABLE_NVECTOR_CUDA OR SUNDIALS_ENABLE_NVECTOR_HIP)
    add_subdirectory(mpi_gpu)
  endif()

  install(FILES README.md
          DESTINATION "${SUNDIALS_BENCHMARKS_INSTALL_PATH}/diffusion_2D")

endif()
