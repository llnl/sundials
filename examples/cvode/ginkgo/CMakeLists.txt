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

# Example lists are tuples "name\;args\;type"

# Examples that support CPU and GPU Ginkgo backends
set(cpu_gpu_examples "cv_heat2D_ginkgo.cpp\;\;develop")

if(SUNDIALS_INDEX_SIZE STREQUAL "32")
  list(APPEND cpu_gpu_examples "cv_bruss_batched_ginkgo.cpp\;\;develop")
endif()

sundials_add_examples_ginkgo(
  cpu_gpu_examples
  TARGETS sundials_cvode
  BACKENDS REF OMP CUDA HIP SYCL)

# Examples that only support CPU Ginkgo backends
set(cpu_examples "cv_kpr_ginkgo.cpp\;\;develop")

sundials_add_examples_ginkgo(
  cpu_examples
  TARGETS sundials_cvode
  BACKENDS REF OMP)

# Install the targets
if(EXAMPLES_INSTALL)

  if(SUNDIALS_GINKGO_BACKENDS MATCHES "CUDA")
    list(APPEND vectors nveccuda)
  endif()
  if(SUNDIALS_GINKGO_BACKENDS MATCHES "HIP")
    list(APPEND vectors nvechip)
  endif()
  if(SUNDIALS_GINKGO_BACKENDS MATCHES "SYCL")
    list(APPEND vectors nvecsycl)
  endif()
  if((SUNDIALS_GINKGO_BACKENDS MATCHES "OMP") OR (SUNDIALS_GINKGO_BACKENDS
                                                  MATCHES "REF"))
    list(APPEND vectors nvecserial)
  endif()

  sundials_install_examples_ginkgo(
    cvode
    CPU_EXAMPLES_VAR cpu_examples
    CPU_GPU_EXAMPLES_VAR cpu_gpu_examples
    SUNDIALS_COMPONENTS cvode ${vectors}
    SUNDIALS_TARGETS cvode
    EXTRA_FILES "${PROJECT_SOURCE_DIR}/examples/utilities/example_utilities.hpp"
    DESTINATION cvode/ginkgo)

endif()
