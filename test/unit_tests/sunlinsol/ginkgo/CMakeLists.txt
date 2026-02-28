# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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

# Example lists are tuples "name\;args\;type" where the type is 'develop' for
# examples excluded from 'make test' in releases
set(examples
    "test_sunlinsol_ginkgo.cpp\;bicg csr 100 0 100 0\;"
    "test_sunlinsol_ginkgo.cpp\;bicg dense 100 0 100 0\;"
    "test_sunlinsol_ginkgo.cpp\;cgs csr 100 0 100 0\;"
    "test_sunlinsol_ginkgo.cpp\;cgs dense 100 0 100 0\;"
    "test_sunlinsol_ginkgo.cpp\;gmres csr 100 0 100 0\;"
    "test_sunlinsol_ginkgo.cpp\;gmres dense 100 0 100 0\;"
    "test_sunlinsol_ginkgo.cpp\;idr csr 100 0 100 0\;"
    "test_sunlinsol_ginkgo.cpp\;idr dense 100 0 100 0\;"
    "test_sunlinsol_ginkgo.cpp\;bicgstab csr 100 0 100 0\;"
    "test_sunlinsol_ginkgo.cpp\;bicgstab dense 100 0 100 0\;"
    "test_sunlinsol_ginkgo_batch.cpp\;bicgstab dense 50 10 2 500 0\;")

if(SUNDIALS_INDEX_SIZE STREQUAL "32")
  list(APPEND examples
       "test_sunlinsol_ginkgo_batch.cpp\;bicgstab csr 50 10 2 500 0\;")
endif()

# Add source directory to include directories
include_directories(..)

sundials_add_examples_ginkgo(
  examples
  TARGETS test_sunlinsol_obj
  BACKENDS REF OMP CUDA HIP SYCL
  UNIT_TEST)
