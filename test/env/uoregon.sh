#!/bin/bash
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
# Environment overrides for the University of Oregon GitLab runners.
# ------------------------------------------------------------------------------

export CMAKE_BUILD_TYPE="${BUILD_TYPE:-RelWithDebInfo}"

# Disable TPLs by default. GPU jobs selectively enable CUDA, HIP, or SYCL in
# the GitLab configuration before this file is sourced.
export SUNDIALS_PTHREAD=OFF
export SUNDIALS_OPENMP=OFF
export SUNDIALS_OPENMP_OFFLOAD=OFF
export SUNDIALS_MPI=OFF
export SUNDIALS_KOKKOS=OFF
export SUNDIALS_RAJA=OFF
export SUNDIALS_GINKGO=OFF
export SUNDIALS_LAPACK=OFF
export SUNDIALS_KLU=OFF
export SUNDIALS_KOKKOS_KERNELS=OFF
export SUNDIALS_SUPERLU_MT=OFF
export SUNDIALS_SUPERLU_DIST=OFF
export SUNDIALS_MAGMA=OFF
export SUNDIALS_HYPRE=OFF
export SUNDIALS_PETSC=OFF
export SUNDIALS_TRILINOS=OFF
export SUNDIALS_XBRAID=OFF
export SUNDIALS_FUSED_KERNELS=OFF

export SUNDIALS_CUDA="${SUNDIALS_CUDA:-OFF}"
export SUNDIALS_HIP="${SUNDIALS_HIP:-OFF}"
export SUNDIALS_SYCL="${SUNDIALS_SYCL:-OFF}"

return 0
