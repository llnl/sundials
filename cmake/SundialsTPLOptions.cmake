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
# SUNDIALS options for third-party libraries
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# Enable MPI support?
# ---------------------------------------------------------------

sundials_option(ENABLE_MPI BOOL "Enable MPI support" OFF DEPRECATED_NAMES
                MPI_ENABLE)

# ---------------------------------------------------------------
# Enable OpenMP support?
# ---------------------------------------------------------------

sundials_option(ENABLE_OPENMP BOOL "Enable OpenMP support" OFF DEPRECATED_NAMES
                OPENMP_ENABLE)

# ---------------------------------------------------------------
# Enable OpenMP target offloading support?
# ---------------------------------------------------------------

sundials_option(
  ENABLE_OPENMP_DEVICE BOOL "Enable OpenMP device offloading support" OFF
  DEPRECATED_NAMES OPENMP_DEVICE_ENABLE)

# Advanced option to skip OpenMP device offloading support check. This is needed
# for a specific compiler that doesn't correctly report its OpenMP spec date
# (with CMake >= 3.9).
sundials_option(OPENMP_DEVICE_WORKS BOOL
                "Skip the OpenMP device offloading support check" OFF ADVANCED)

# ---------------------------------------------------------------
# Enable Pthread support?
# ---------------------------------------------------------------

sundials_option(ENABLE_PTHREAD BOOL "Enable Pthreads support" OFF
                DEPRECATED_NAMES PTHREAD_ENABLE)

# -------------------------------------------------------------
# Enable CUDA support?
# -------------------------------------------------------------

sundials_option(ENABLE_CUDA BOOL "Enable CUDA support" OFF DEPRECATED_NAMES
                CUDA_ENABLE)

# -------------------------------------------------------------
# Enable HIP support?
# -------------------------------------------------------------

sundials_option(ENABLE_HIP BOOL "Enable HIP support" OFF)

# -------------------------------------------------------------
# Enable SYCL support?
# -------------------------------------------------------------

sundials_option(ENABLE_SYCL BOOL "Enable SYCL support" OFF)

sundials_option(
  SUNDIALS_SYCL_2020_UNSUPPORTED
  BOOL
  "Disable the use of some SYCL 2020 features in SUNDIALS libraries and examples"
  OFF
  DEPENDS_ON ENABLE_SYCL
  ADVANCED)

# ---------------------------------------------------------------
# Enable LAPACK support?
# ---------------------------------------------------------------

sundials_option(ENABLE_LAPACK BOOL "Enable Lapack support" OFF DEPRECATED
                LAPACK_ENABLE)

sundials_option(LAPACK_LIBRARIES STRING "Lapack and Blas libraries"
                "${LAPACK_LIBRARIES}")

sundials_option(
  LAPACK_WORKS BOOL
  "Set to ON to force CMake to accept a given LAPACK configuration" OFF
  ADVANCED)

# ---------------------------------------------------------------
# Enable Ginkgo support?
# ---------------------------------------------------------------

sundials_option(ENABLE_GINKGO BOOL "Enable Ginkgo support" OFF)

sundials_option(Ginkgo_DIR PATH "Path to the root of a Ginkgo installation"
                "${Ginkgo_DIR}")

sundials_option(
  SUNDIALS_GINKGO_BACKENDS
  STRING
  "Which Ginkgo backend(s) to build the SUNDIALS Ginkgo interfaces for (REF, OMP, CUDA, HIP, SYCL)"
  "REF;OMP"
  DEPENDS_ON ENABLE_GINKGO)

sundials_option(
  GINKGO_WORKS BOOL
  "Set to ON to force CMake to accept a given Ginkgo configuration" OFF
  ADVANCED)

# ---------------------------------------------------------------
# Enable MAGMA support?
# ---------------------------------------------------------------

sundials_option(ENABLE_MAGMA BOOL "Enable MAGMA support" OFF)

sundials_option(MAGMA_DIR PATH "Path to the root of a MAGMA installation"
                "${MAGMA_DIR}")

sundials_option(
  SUNDIALS_MAGMA_BACKENDS STRING
  "Which MAGMA backend to use under the SUNDIALS MAGMA interfaces (CUDA, HIP)"
  "CUDA"
  OPTIONS "CUDA;HIP"
  DEPENDS_ON ENABLE_MAGMA)

sundials_option(
  MAGMA_WORKS BOOL
  "Set to ON to force CMake to accept a given MAGMA configuration" OFF ADVANCED)

# ---------------------------------------------------------------
# Enable SuperLU_DIST support?
# ---------------------------------------------------------------

sundials_option(ENABLE_SUPERLUDIST BOOL "Enable SuperLU_DIST support" OFF
                DEPRECATED_NAMES SUPERLUDIST_ENABLE)

sundials_option(
  SUPERLUDIST_DIR PATH "Path to the root of the SuperLU_DIST installation"
  "${SUPERLUDIST_DIR}")

sundials_option(
  SUPERLUDIST_INCLUDE_DIRS
  PATH
  "SuperLU_DIST include directories"
  "${SUPERLUDIST_INCLUDE_DIRS}"
  ADVANCED
  DEPRECATED_NAMES
  SUPERLUDIST_INCLUDE_DIR)

sundials_option(
  SUPERLUDIST_LIBRARIES STRING
  "Semi-colon separated list of libraries needed for SuperLU_DIST."
  "${SUPERLUDIST_LIBRARIES}" ADVANCED)

sundials_option(
  SUPERLUDIST_OpenMP BOOL
  "Enable SUNDIALS support for SuperLU_DIST OpenMP on-node parallelism" OFF)

sundials_option(
  SUPERLUDIST_WORKS BOOL
  "Set to ON to force CMake to accept a given SuperLU_DIST configuration" OFF
  ADVANCED)

# ---------------------------------------------------------------
# Enable SuperLU_MT support?
# ---------------------------------------------------------------

sundials_option(ENABLE_SUPERLUMT BOOL "Enable SuperLU_MT support" OFF
                DEPRECATED_NAMES SUPERLUMT_ENABLE)

sundials_option(SUPERLUMT_INCLUDE_DIR PATH "SuperLU_MT include directory"
                "${SUPERLUMT_INCLUDE_DIR}")

sundials_option(SUPERLUMT_LIBRARY_DIR PATH "SuperLU_MT library directory"
                "${SUPERLUMT_LIBRARY_DIR}")

sundials_option(
  SUPERLUMT_LIBRARIES STRING
  "Semi-colon separated list of additional libraries needed for SuperLU_MT."
  "${SUPERLUMT_LIBRARIES}")

sundials_option(SUPERLUMT_THREAD_TYPE STRING
                "SuperLU_MT threading type: OPENMP or PTHREAD" "PTHREAD")

sundials_option(
  SUPERLUMT_WORKS BOOL
  "Set to ON to force CMake to accept a given SUPERLUMT configuration" OFF
  ADVANCED)

# ---------------------------------------------------------------
# Enable KLU support?
# ---------------------------------------------------------------

sundials_option(ENABLE_KLU BOOL "Enable KLU support" OFF DEPRECATED_NAMES
                KLU_ENABLE)

sundials_option(KLU_INCLUDE_DIR PATH "KLU include directory"
                "${KLU_INCLUDE_DIR}")

sundials_option(KLU_LIBRARY_DIR PATH "KLU library directory"
                "${KLU_LIBRARY_DIR}")

sundials_option(
  KLU_WORKS BOOL "Set to ON to force CMake to accept a given KLU configuration"
  OFF ADVANCED)

# ---------------------------------------------------------------
# Enable hypre support?
# ---------------------------------------------------------------

sundials_option(ENABLE_HYPRE BOOL "Enable hypre support" OFF DEPRECATED_NAMES
                HYPRE_ENABLE)

sundials_option(HYPRE_DIR PATH "Path to hypre installation" "${HYPRE_DIR}")

sundials_option(HYPRE_INCLUDE_DIR PATH "HYPRE include directory"
                "${HYPRE_INCLUDE_DIR}")

sundials_option(HYPRE_LIBRARY_DIR PATH "HYPRE library directory"
                "${HYPRE_LIBRARY_DIR}")

sundials_option(
  HYPRE_WORKS BOOL
  "Set to ON to force CMake to accept a given hypre configuration" OFF ADVANCED)

# ---------------------------------------------------------------
# Enable PETSc support?
# ---------------------------------------------------------------

sundials_option(ENABLE_PETSC BOOL "Enable PETSc support" OFF DEPRECATED_NAMES
                PETSC_ENABLE)

sundials_option(PETSC_DIR PATH "Path to the root of a PETSc installation"
                "${PETSC_DIR}")

sundials_option(PETSC_ARCH STRING "PETSc architecture (optional)"
                "${PETSC_ARCH}")

sundials_option(
  PETSC_LIBRARIES STRING "Semi-colon separated list of PETSc link libraries"
  "${PETSC_LIBRARIES}" ADVANCED)

sundials_option(
  PETSC_INCLUDES STRING
  "Semi-colon separated list of PETSc include directories" "${PETSC_INCLUDES}"
  ADVANCED)

sundials_option(
  PETSC_WORKS BOOL
  "Set to ON to force CMake to accept a given PETSc configuration" OFF ADVANCED)

# -------------------------------------------------------------
# Enable RAJA support?
# -------------------------------------------------------------

sundials_option(ENABLE_RAJA BOOL "Enable RAJA support" OFF DEPRECATED_NAMES
                RAJA_ENABLE)

sundials_option(RAJA_DIR PATH "Path to root of RAJA installation" "${RAJA_DIR}")

sundials_option(
  SUNDIALS_RAJA_BACKENDS STRING
  "Which RAJA backend under the SUNDIALS RAJA interfaces (CUDA, HIP, SYCL)"
  "CUDA"
  OPTIONS "CUDA;HIP;SYCL"
  DEPENDS_ON ENABLE_RAJA)

# ---------------------------------------------------------------
# Enable Trilinos support?
# ---------------------------------------------------------------

sundials_option(ENABLE_TRILINOS BOOL "Enable Trilinos support" OFF
                DEPRECATED_NAMES Trilinos_ENABLE)

sundials_option(Trilinos_DIR PATH "Path to root of Trilinos installation"
                "${Trilinos_DIR}")

# ---------------------------------------------------------------
# Enable XBraid support?
# ---------------------------------------------------------------

sundials_option(ENABLE_XBRAID BOOL "Enable XBraid support" OFF)

sundials_option(XBRAID_DIR PATH "Path to the root of an XBraid installation"
                "${XBRAID_DIR}")

sundials_option(
  XBRAID_LIBRARIES STRING "Semi-colon separated list of XBraid link libraries"
  "${XBRAID_LIBRARIES}" ADVANCED)

sundials_option(
  XBRAID_INCLUDES STRING
  "Semi-colon separated list of XBraid include directories"
  "${XBRAID_INCLUDES}" ADVANCED)

sundials_option(
  XBRAID_WORKS BOOL
  "Set to ON to force CMake to accept a given XBraid configuration" OFF
  ADVANCED)

# -------------------------------------------------------------
# Enable oneMKL support?
# -------------------------------------------------------------

sundials_option(ENABLE_ONEMKL BOOL "Enable oneMKL support" OFF)

sundials_option(ONEMKL_DIR PATH "Path to root of oneMKL installation"
                "${ONEMKL_DIR}")

sundials_option(
  ONEMKL_WORKS BOOL
  "Set to ON to force CMake to accept a given oneMKL configuration" OFF
  ADVANCED)

sundials_option(
  SUNDIALS_ONEMKL_USE_GETRF_LOOP BOOL
  "Replace batched getrf call with loop over getrf" OFF
  DEPENDS_ON ENABLE_ONEMKL
  ADVANCED)

sundials_option(
  SUNDIALS_ONEMKL_USE_GETRS_LOOP BOOL
  "Replace batched getrs call with loop over getrs" OFF
  DEPENDS_ON ENABLE_ONEMKL
  ADVANCED)

# ---------------------------------------------------------------
# Enable Caliper support?
# ---------------------------------------------------------------

sundials_option(ENABLE_CALIPER BOOL "Enable CALIPER support" OFF
                DEPENDS_ON SUNDIALS_BUILD_WITH_PROFILING)

sundials_option(CALIPER_DIR PATH "Path to the root of an CALIPER installation"
                "${CALIPER_DIR}")

sundials_option(
  CALIPER_WORKS BOOL
  "Set to ON to force CMake to accept a given CALIPER configuration" OFF
  ADVANCED)

# ---------------------------------------------------------------
# Enable Adiak support?
# ---------------------------------------------------------------

sundials_option(ENABLE_ADIAK BOOL "Enable Adiak support" OFF
                DEPENDS_ON SUNDIALS_BUILD_WITH_PROFILING)

sundials_option(adiak_DIR PATH "Path to the root of an Adiak installation"
                "${ADIAK_DIR}")

# ---------------------------------------------------------------
# Enable Kokkos support?
# ---------------------------------------------------------------

sundials_option(ENABLE_KOKKOS BOOL "Enable Kokkos support" OFF)

sundials_option(Kokkos_DIR PATH "Path to the root of a Kokkos installation"
                "${Kokkos_DIR}")

sundials_option(
  KOKKOS_WORKS BOOL
  "Set to ON to force CMake to accept a given Kokkos configuration" OFF
  ADVANCED)

# ---------------------------------------------------------------
# Enable Kokkos Kernels support?
# ---------------------------------------------------------------

sundials_option(ENABLE_KOKKOS_KERNELS BOOL "Enable Kokkos Kernels support" OFF)

sundials_option(
  KokkosKernels_DIR PATH "Path to the root of a Kokkos Kernels installation"
  "${KokkosKernels_DIR}")

sundials_option(
  KOKKOS_KERNELS_WORKS BOOL
  "Set to ON to force CMake to accept a given Kokkos configuration" OFF
  ADVANCED)
