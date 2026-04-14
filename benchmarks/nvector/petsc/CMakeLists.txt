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
# CMakeLists.txt file for PETSc nvector examples
# ---------------------------------------------------------------

message(STATUS "Added PETSc NVECTOR benchmark")

sundials_add_nvector_benchmark(
  nvector_petsc_benchmark
  SOURCES test_nvector_performance_petsc.c
  SUNDIALS_TARGETS sundials_nvecpetsc
  LINK_LIBRARIES MPI::MPI_C PUBLIC SUNDIALS::PETSC
  INSTALL_SUBDIR nvector/petsc)
