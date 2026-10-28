# -----------------------------------------------------------------
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
# -----------------------------------------------------------------
# CVODE Ginkgo batched Brusselator example with CPU backends
set(cpu_examples "cv_bruss_batched_ginkgo.cpp\;\;develop")
sundials_add_examples_ginkgo(
  cpu_examples
  TARGETS ${CVODE_C_LIBS}
  BACKENDS REF OMP INSTALL_FILES cv_bruss_batched_ginkgo.hpp
           "${PROJECT_SOURCE_DIR}/examples/utilities/example_utilities.hpp")
