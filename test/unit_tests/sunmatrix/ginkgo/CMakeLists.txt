# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# ---------------------------------------------------------------
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
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is 'develop' for
# examples excluded from 'make test' in releases
set(examples
    "test_sunmatrix_ginkgo.cpp\;100 100 0\;"
    "test_sunmatrix_ginkgo.cpp\;100 10 0\;"
    "test_sunmatrix_ginkgo.cpp\;10 100 0\;"
    "test_sunmatrix_ginkgo.cpp\;100 100 1\;"
    "test_sunmatrix_ginkgo.cpp\;100 10 1\;"
    "test_sunmatrix_ginkgo.cpp\;10 100 1\;"
    "test_sunmatrix_ginkgo_batch.cpp\;10 10 10 1\;"
    "test_sunmatrix_ginkgo_batch.cpp\;10 100 10 1\;")

if(SUNDIALS_INDEX_SIZE STREQUAL "32")
  list(APPEND examples "test_sunmatrix_ginkgo_batch.cpp\;10 10 10 0\;"
       "test_sunmatrix_ginkgo_batch.cpp\;10 100 10 0\;")
endif()

# Add source directory to include directories
include_directories(..)

sundials_add_examples_ginkgo(
  examples
  TARGETS test_sunmatrix_obj sundials_sunmatrixdense
  BACKENDS REF OMP CUDA HIP SYCL
  UNIT_TEST)
