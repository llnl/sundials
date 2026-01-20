#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# -----------------------------------------------------------------------------
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
# -----------------------------------------------------------------------------
# Script to setup a CMake cache file or configure script for SUNDIALS
# -----------------------------------------------------------------------------


def main():
    """Create a CMake cache file or configure script"""

    import argparse

    parser = argparse.ArgumentParser(description="""Create a SUNDIALS CMake
                                     cache file""")

    parser.add_argument(
        "--filetype",
        type=str,
        choices=["cache", "script"],
        default="cache",
        help="""Create a CMake cache file or configuration
                        script (default cache)""",
    )

    parser.add_argument(
        "--filename",
        type=str,
        default="sundials.cmake",
        help="""Set the cache file or script name (default
                        sundials.cmake)""",
    )

    parser.add_argument(
        "--readenv",
        action="store_true",
        help="""Read environment variables (command line
                        arguments will override any settings from the
                        environment variables)""",
    )

    parser.add_argument(
        "--debugscript", action="store_true", help="Enable debugging output for this script"
    )

    # -----------------
    # Compiler Options
    # -----------------

    group = parser.add_argument_group(
        "Compilers and Flags",
        """Options for setting the C, C++,
                                      Fortran, and CUDA compiler and flags.""",
    )

    # Build type
    add_arg(
        group,
        "--build-type",
        "CMAKE_BUILD_TYPE",
        "CMAKE_BUILD_TYPE",
        "RelWithDebInfo",
        "STRING",
        "CMake build type (Debug, RelWithDebInfo, Release)",
    )

    # C compiler
    add_arg(group, "--c-compiler", "CC", "CMAKE_C_COMPILER", None, "FILEPATH", "C compiler")

    add_arg(group, "--c-flags", "CFLAGS", "CMAKE_C_FLAGS", None, "STRING", "C compiler flags")

    add_arg(group, "--c-std", "CMAKE_C_STANDARD", "CMAKE_C_STANDARD", "99", "STRING", "C standard")

    add_arg(
        group,
        "--c-ext",
        "CMAKE_C_EXTENSIONS",
        "CMAKE_C_EXTENSIONS",
        "OFF",
        "STRING",
        "C compiler extensions",
    )

    # C++ compiler
    add_arg(group, "--cxx-compiler", "CXX", "CMAKE_CXX_COMPILER", None, "FILEPATH", "C++ compiler")

    add_arg(
        group, "--cxx-flags", "CXXFLAGS", "CMAKE_CXX_FLAGS", None, "STRING", "C++ compiler flags"
    )

    add_arg(
        group,
        "--cxx-std",
        "CMAKE_CXX_STANDARD",
        "CMAKE_CXX_STANDARD",
        "14",
        "STRING",
        "C++ standard",
    )

    add_arg(
        group,
        "--cxx-ext",
        "CMAKE_CXX_EXTENSIONS",
        "CMAKE_CXX_EXTENSIONS",
        "OFF",
        "STRING",
        "C++ compiler extensions",
    )

    # Fortran compiler
    add_arg(
        group,
        "--fortran-compiler",
        "FC",
        "CMAKE_Fortran_COMPILER",
        None,
        "FILEPATH",
        "Fortran compiler",
    )

    add_arg(
        group,
        "--fortran-flags",
        "FFLAGS",
        "CMAKE_Fortran_FLAGS",
        None,
        "STRING",
        "Fortran compiler flags",
    )

    # CUDA compiler
    add_arg(
        group,
        "--cuda-compiler",
        "CUDACXX",
        "CMAKE_CUDA_COMPILER",
        None,
        "FILEPATH",
        "CUDA compiler",
    )

    add_arg(
        group,
        "--cuda-flags",
        "CUDAFLAGS",
        "CMAKE_CUDA_FLAGS",
        None,
        "STRING",
        "CUDA compiler flags",
    )

    add_arg(
        group,
        "--cuda-std",
        "CMAKE_CUDA_STANDARD",
        "CMAKE_CUDA_STANDARD",
        "14",
        "STRING",
        "CUDA standard",
    )

    add_arg(
        group,
        "--cuda-arch",
        "CUDAARCHS",
        "CMAKE_CUDA_ARCHITECTURES",
        None,
        "STRING",
        "CUDA architecture",
    )

    # Additional compiler options
    add_arg(
        group,
        "--Wall",
        "SUNDIALS_ENABLE_ALL_WARNINGS",
        "ENABLE_ALL_WARNINGS",
        "OFF",
        "BOOL",
        "Enable all compiler warnings",
    )

    add_arg(
        group,
        "--Werror",
        "SUNDIALS_ENABLE_WARNINGS_AS_ERRORS",
        "ENABLE_WARNINGS_AS_ERRORS",
        "OFF",
        "BOOL",
        "Enable compiler warnings as errors",
    )

    add_arg(
        group,
        "--address-sanitizer",
        "SUNDIALS_ENABLE_ADDRESS_SANITIZER",
        "ENABLE_ADDRESS_SANITIZER",
        "OFF",
        "BOOL",
        "Enable address sanitizer",
    )

    # ----------------
    # Install Options
    # ----------------

    group = parser.add_argument_group(
        "Install Options",
        """Options for where SUNDIALS should be
                                      installed.""",
    )

    # install prefix
    add_arg(
        group,
        "--install-prefix",
        "SUNDIALS_INSTALL_PREFIX",
        "CMAKE_INSTALL_PREFIX",
        None,
        "PATH",
        "SUNDIALS install location",
    )

    # library directory

    # ------------------
    # Debugging Options
    # ------------------

    group = parser.add_argument_group("Debugging Options", """Options debugging SUNDIALS.""")

    add_arg(
        group,
        "--debug",
        "SUNDIALS_DEBUG",
        "SUNDIALS_DEBUG",
        "OFF",
        "BOOL",
        "SUNDIALS debugging output",
    )

    add_arg(
        group,
        "--debug-assert",
        "SUNDIALS_DEBUG_ASSERT",
        "SUNDIALS_DEBUG_ASSERT",
        "OFF",
        "BOOL",
        "SUNDIALS debugging asserts",
        dependson="--debug",
    )

    add_arg(
        group,
        "--debug-cuda",
        "SUNDIALS_DEBUG_CUDA_LASTERROR",
        "SUNDIALS_DEBUG_CUDA_LASTERROR",
        "OFF",
        "BOOL",
        "SUNDIALS debugging cuda errors",
        dependson="--debug",
    )

    add_arg(
        group,
        "--debug-hip",
        "SUNDIALS_DEBUG_HIP_LASTERROR",
        "SUNDIALS_DEBUG_HIP_LASTERROR",
        "OFF",
        "BOOL",
        "SUNDIALS debugging hip errors",
        dependson="--debug",
    )

    add_arg(
        group,
        "--debug-printvec",
        "SUNDIALS_DEBUG_PRINTVEC",
        "SUNDIALS_DEBUG_PRINTVEC",
        "OFF",
        "BOOL",
        "SUNDIALS debugging vector output",
        dependson="--debug",
    )

    # --------------
    # Library Types
    # --------------

    group = parser.add_argument_group(
        "Library Type Options",
        """Options to specify if shared and/or
                                      static libraries are build.""",
    )
    add_arg(
        group,
        "--static",
        "SUNDIALS_STATIC_LIBRARIES",
        "BUILD_STATIC_LIBS",
        "ON",
        "BOOL",
        "Build static SUNDIALS libraries",
    )

    add_arg(
        group,
        "--shared",
        "SUNDIALS_SHARED_LIBRARIES",
        "BUILD_SHARED_LIBS",
        "ON",
        "BOOL",
        "Build shared SUNDIALS libraries",
    )

    # ---------
    # Packages
    # ---------

    # packages TODO(DJG): Add support for ONLY option
    group = parser.add_argument_group(
        "SUNDIALS Packages",
        """Options to specify which SUNDIALS
                                      packages should be built.""",
    )
    add_arg(
        group,
        "--arkode",
        "SUNDIALS_ARKODE",
        "BUILD_ARKODE",
        "ON",
        "BOOL",
        "Build the ARKODE library",
    )

    add_arg(
        group, "--cvode", "SUNDIALS_CVODE", "BUILD_CVODE", "ON", "BOOL", "Build the CVODE library"
    )

    add_arg(
        group,
        "--cvodes",
        "SUNDIALS_CVODES",
        "BUILD_CVODES",
        "ON",
        "BOOL",
        "Build the CVODES library",
    )

    add_arg(group, "--ida", "SUNDIALS_IDA", "BUILD_IDA", "ON", "BOOL", "Build the IDA library")

    add_arg(group, "--idas", "SUNDIALS_IDAS", "BUILD_IDAS", "ON", "BOOL", "Build the IDAS library")

    add_arg(
        group,
        "--kinsol",
        "SUNDIALS_KINSOL",
        "BUILD_KINSOL",
        "ON",
        "BOOL",
        "Build the KINSOL library",
    )

    # -----------------
    # Packages Options
    # -----------------

    group = parser.add_argument_group(
        "SUNDIALS Package Options",
        """Options for configuring SUNDIALS types
                                      and enabling special compile time
                                      features.""",
    )
    # index size
    add_arg(
        group,
        "--indexsize",
        "SUNDIALS_INDEX_SIZE",
        "SUNDIALS_INDEX_SIZE",
        "64",
        "STRING",
        "index size",
        choices=["32", "64"],
    )

    # precision
    add_arg(
        group,
        "--precision",
        "SUNDIALS_PRECISION",
        "SUNDIALS_PRECISION",
        "double",
        "STRING",
        "real type precision",
        choices=["single", "double", "extended"],
    )

    # monitoring
    add_arg(
        group,
        "--monitoring",
        "SUNDIALS_MONITORING",
        "SUNDIALS_BUILD_WITH_MONITORING",
        "OFF",
        "BOOL",
        "integrator and solver monitoring",
    )

    # profiling
    add_arg(
        group,
        "--profiling",
        "SUNDIALS_PROFILING",
        "SUNDIALS_BUILD_WITH_PROFILING",
        "OFF",
        "BOOL",
        "fine-grained profiling",
    )

    add_arg(
        group,
        "--logging-level",
        "SUNDIALS_LOGGING_LEVEL",
        "SUNDIALS_LOGGING_LEVEL",
        "0",
        "STRING",
        "logging",
        choices=["0", "1", "2", "3", "4", "5"],
    )

    # fused kernels
    add_arg(
        group,
        "--fused-kernels",
        "SUNDIALS_FUSED_KERNELS",
        "SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS",
        "OFF",
        "BOOL",
        "package fused kernels",
    )

    # error checks
    add_arg(
        group,
        "--enable-error-checks",
        "SUNDIALS_ENABLE_ERROR_CHECKS",
        "SUNDIALS_ENABLE_ERROR_CHECKS",
        "OFF",
        "BOOL",
        "enable error checks",
    )

    # -----------
    # Interfaces
    # -----------

    group = parser.add_argument_group(
        "SUNDIALS Interfaces",
        """These options enable or disable the
                                      SUNDIALS Fortran interfaces.""",
    )

    # Fortran interfaces
    add_arg(
        group,
        "--fmod-interface",
        "SUNDIALS_FMOD_INTERFACE",
        "BUILD_FORTRAN_MODULE_INTERFACE",
        "OFF",
        "BOOL",
        "Fortran module interface",
    )

    # ---------
    # Examples
    # ---------

    group = parser.add_argument_group(
        "Example and Benchmark Programs",
        """These options enable or disable
                                      building and installing the SUNDIALS
                                      example and Benchmark programs.""",
    )

    add_arg(
        group,
        "--examples-c",
        "SUNDIALS_EXAMPLES_C",
        "EXAMPLES_ENABLE_C",
        "ON",
        "BOOL",
        "C examples",
    )

    add_arg(
        group,
        "--examples-cxx",
        "SUNDIALS_EXAMPLES_CXX",
        "EXAMPLES_ENABLE_CXX",
        None,
        "BOOL",
        "C++ examples",
    )

    add_arg(
        group,
        "--examples-f03",
        "SUNDIALS_EXAMPLES_F03",
        "EXAMPLES_ENABLE_F2003",
        None,
        "BOOL",
        "Fortran 2003 examples",
    )

    add_arg(
        group,
        "--examples-cuda",
        "SUNDIALS_EXAMPLES_CUDA",
        "EXAMPLES_ENABLE_CUDA",
        None,
        "BOOL",
        "CUDA examples",
    )

    add_arg(
        group,
        "--benchmarks",
        "SUNDIALS_BENCHMARKS",
        "BUILD_BENCHMARKS",
        "OFF",
        "BOOL",
        "Benchmarks",
    )

    # ------------
    # TPL Options
    # ------------

    # ----
    # MPI
    # ----

    group = parser.add_argument_group(
        "MPI Options",
        """Options for enabling MPI support in
                                      SUNDIALS and setting the MPI C, C++, and
                                      Fortran compilers.""",
    )

    add_arg(
        group, "--mpi", "SUNDIALS_MPI", "ENABLE_MPI", "OFF", "FILEPATH", "SUNDIALS MPI support"
    )

    add_arg(
        group,
        "--mpicc",
        "MPICC",
        "MPI_C_COMPILER",
        None,
        "FILEPATH",
        "MPI C compiler",
        dependson="--mpi",
    )

    add_arg(
        group,
        "--mpicxx",
        "MPICXX",
        "MPI_CXX_COMPILER",
        None,
        "FILEPATH",
        "MPI C++ compiler",
        dependson="--mpi",
    )

    add_arg(
        group,
        "--mpifort",
        "MPIFC",
        "MPI_Fortran_COMPILER",
        None,
        "FILEPATH",
        "MPI Fortran compiler",
        dependson="--mpi",
    )

    add_arg(
        group,
        "--mpiexec",
        "MPIEXEC",
        "MPIEXEC_EXECUTABLE",
        None,
        "FILEPATH",
        "MPI executable",
        dependson="--mpi",
    )

    add_arg(
        group,
        "--mpiexec-pre-flags",
        "MPIEXEC_PREFLAGS",
        "MPIEXEC_PREFLAGS",
        None,
        "STRING",
        "MPI executable extra flags",
        dependson="--mpi",
    )

    # ----------
    # Threading
    # ----------

    # OpenMP
    group = parser.add_argument_group(
        "OpenMP Options",
        """Options for enabling OpenMP support in
                                      SUNDIALS.""",
    )

    add_arg(
        group,
        "--openmp",
        "SUNDIALS_OPENMP",
        "ENABLE_OPENMP",
        "OFF",
        "BOOL",
        "SUNDIALS OpenMP support",
    )

    add_arg(
        group,
        "--openmp-device-works",
        "SUNDIALS_OPENMP_DEVICE_WORKS",
        "OPENMP_DEVICE_WORKS",
        "OFF",
        "BOOL",
        "Disable OpenMP Device Support Checks (assume OpenMP 4.5+)",
    )

    # Pthread
    group = parser.add_argument_group(
        "Pthread Options",
        """Options for enabling
                                      Pthread support in SUNDIALS.""",
    )

    add_arg(
        group,
        "--pthread",
        "SUNDIALS_PTHREAD",
        "ENABLE_PTHREAD",
        "OFF",
        "BOOL",
        "SUNDIALS PThread support",
    )

    # -----
    # GPUs
    # -----

    # CUDA
    group = parser.add_argument_group(
        "CUDA Options",
        """Options for enabling CUDA support in

                                      SUNDIALS""",
    )
    add_arg(
        group, "--cuda", "SUNDIALS_CUDA", "ENABLE_CUDA", "OFF", "BOOL", "SUNDIALS CUDA support"
    )

    # HIP
    group = parser.add_argument_group(
        "HIP Options",
        """Options for enabling HIP support in
                                      SUNDIALS.""",
    )

    add_arg(group, "--hip", "SUNDIALS_HIP", "ENABLE_HIP", "OFF", "BOOL", "SUNDIALS HIP support")

    # OpenMP Offload
    group = parser.add_argument_group(
        "OpenMP Offload Options",
        """Options for enabling OpenMP offload
                                      support in SUNDIALS.""",
    )

    add_arg(
        group,
        "--openmp-offload",
        "SUNDIALS_OPENMP_OFFLOAD",
        "ENABLE_OPENMP_DEVICE",
        "OFF",
        "BOOL",
        "SUNDIALS OpenMP offload support",
    )

    # ------------------------
    # Performance portability
    # ------------------------

    # Kokkos
    group = parser.add_argument_group("Kokkos Options")

    add_arg(
        group,
        "--kokkos",
        "SUNDIALS_KOKKOS",
        "ENABLE_KOKKOS",
        "OFF",
        "BOOL",
        "SUNDIALS Kokkos support",
    )

    add_arg(
        group,
        "--kokkos-dir",
        "KOKKOS_ROOT",
        "Kokkos_DIR",
        None,
        "PATH",
        "Kokkos install directory",
        dependson="--kokkos",
    )

    # RAJA
    group = parser.add_argument_group("RAJA Options")

    add_arg(
        group, "--raja", "SUNDIALS_RAJA", "ENABLE_RAJA", "OFF", "BOOL", "SUNDIALS Raja support"
    )

    add_arg(
        group,
        "--raja-dir",
        "RAJA_ROOT",
        "RAJA_DIR",
        None,
        "PATH",
        "RAJA install directory",
        dependson="--raja",
    )

    add_arg(
        group,
        "--raja-backends",
        "RAJA_BACKENDS",
        "SUNDIALS_RAJA_BACKENDS",
        None,
        "STRING",
        "RAJA backends",
        choices=["CUDA", "HIP"],
        dependson="--raja",
    )

    # SYCL
    group = parser.add_argument_group("SYCL Options")

    add_arg(
        group, "--sycl", "SUNDIALS_SYCL", "ENABLE_SYCL", "OFF", "BOOL", "SUNDIALS SYCL support"
    )

    # ------------------------
    # Linear solver libraries
    # ------------------------

    # Ginkgo
    group = parser.add_argument_group("Ginkgo Options")

    add_arg(
        group,
        "--ginkgo",
        "SUNDIALS_GINKGO",
        "ENABLE_GINKGO",
        "OFF",
        "BOOL",
        "SUNDIALS Ginkgo support",
    )

    add_arg(
        group,
        "--ginkgo-dir",
        "GINKGO_ROOT",
        "Ginkgo_DIR",
        None,
        "PATH",
        "Ginkgo install directory",
        dependson="--ginkgo",
    )

    add_arg(
        group,
        "--ginkgo-backends",
        "GINKGO_BACKENDS",
        "SUNDIALS_GINKGO_BACKENDS",
        "REF;OMP",
        "STRING",
        "Ginkgo backends",
        choices=["REF", "OMP", "CUDA", "HIP", "DPCPP"],
        dependson="--ginkgo",
    )

    # LAPACK
    group = parser.add_argument_group("LAPACK Options")

    add_arg(
        group,
        "--lapack",
        "SUNDIALS_LAPACK",
        "ENABLE_LAPACK",
        "OFF",
        "BOOL",
        "SUNDIALS LAPACK support",
    )

    add_arg(
        group,
        "--lapack-libs",
        "LAPACK_LIBRARIES",
        "LAPACK_LIBRARIES",
        None,
        "STRING",
        "LAPACK libraries",
        dependson="--lapack",
    )

    add_arg(
        group,
        "--blas-libs",
        "BLAS_LIBRARIES",
        "BLAS_LIBRARIES",
        None,
        "STRING",
        "BLAS libraries",
        dependson="--lapack",
    )

    # KLU
    group = parser.add_argument_group("KLU Options")

    add_arg(group, "--klu", "SUNDIALS_KLU", "ENABLE_KLU", "OFF", "BOOL", "SUNDIALS KLU support")

    add_arg(
        group,
        "--klu-incdir",
        "SUITE_SPARSE_INCLUDE_DIR",
        "KLU_INCLUDE_DIR",
        None,
        "PATH",
        "KLU include directory",
        dependson="--klu",
    )

    add_arg(
        group,
        "--klu-libdir",
        "SUITE_SPARSE_LIBRARY_DIR",
        "KLU_LIBRARY_DIR",
        None,
        "PATH",
        "KLU library directory",
        dependson="--klu",
    )

    # KokkosKernels
    group = parser.add_argument_group("KokkosKernels Options")

    add_arg(
        group,
        "--kokkos-kernels",
        "SUNDIALS_KOKKOS_KERNELS",
        "ENABLE_KOKKOS_KERNELS",
        "OFF",
        "BOOL",
        "SUNDIALS Kokkos-Kernels support",
    )

    add_arg(
        group,
        "--kokkos-kernels-dir",
        "KOKKOS_KERNELS_ROOT",
        "KokkosKernels_DIR",
        None,
        "PATH",
        "Kokkos-Kernels install directory",
        dependson="--kokkos-kernels",
    )

    # SuperLU MT
    group = parser.add_argument_group("SuperLU_MT Options")

    add_arg(
        group,
        "--superlu-mt",
        "SUNDIALS_SUPERLU_MT",
        "ENABLE_SUPERLUMT",
        "OFF",
        "BOOL",
        "SUNDIALS SuperLU MT support",
    )

    add_arg(
        group,
        "--superlu-mt-incdir",
        "SUPERLU_MT_INCLUDE_DIR",
        "SUPERLUMT_INCLUDE_DIR",
        None,
        "PATH",
        "SuperLU_MT include directory",
        dependson="--superlu-mt",
    )

    add_arg(
        group,
        "--superlu-mt-libdir",
        "SUPERLU_MT_LIBRARY_DIR",
        "SUPERLUMT_LIBRARY_DIR",
        None,
        "PATH",
        "SuperLU_MT library directory",
        dependson="--superlu-mt",
    )

    add_arg(
        group,
        "--superlu-mt-libs",
        "SUPERLU_MT_LIBRARIES",
        "SUPERLUMT_LIBRARIES",
        None,
        "STRING",
        "SuperLU_MT additional libraries",
        dependson="--superlu-mt",
    )

    add_arg(
        group,
        "--superlu-mt-thread-type",
        "SUPERLU_MT_THREAD_TYPE",
        "SUPERLUMT_THREAD_TYPE",
        None,
        "STRING",
        "SuperLU_MT thread type",
        choices=["OpenMP", "Pthread"],
        dependson="--superlu-mt",
    )

    # SuperLU DIST
    group = parser.add_argument_group("SuperLU_DIST Options")

    add_arg(
        group,
        "--superlu-dist",
        "SUNDIALS_SUPERLU_DIST",
        "ENABLE_SUPERLUDIST",
        "OFF",
        "BOOL",
        "SUNDIALS SuperLU DIST support",
    )

    add_arg(
        group,
        "--superlu-dist-dir",
        "SUPERLU_DIST_ROOT",
        "SUPERLUDIST_DIR",
        None,
        "PATH",
        "SuperLU_DIST installation directory",
        dependson="--superlu-dist",
    )

    add_arg(
        group,
        "--superlu-dist-incdir",
        "SUPERLU_DIST_INCLUDE_DIR",
        "SUPERLUDIST_INCLUDE_DIR",
        None,
        "PATH",
        "SuperLU_DIST include directory",
        dependson="--superlu-dist",
    )

    add_arg(
        group,
        "--superlu-dist-libdir",
        "SUPERLU_DIST_LIBRARY_DIR",
        "SUPERLUDIST_LIBRARY_DIR",
        None,
        "PATH",
        "SuperLU_DIST library directory",
        dependson="--superlu-dist",
    )

    add_arg(
        group,
        "--superlu-dist-libs",
        "SUPERLU_DIST_LIBRARIES",
        "SUPERLUDIST_LIBRARIES",
        None,
        "STRING",
        "SuperLU_DIST additional libraries",
        dependson="--superlu-dist",
    )

    add_arg(
        group,
        "--superlu-dist-openmp",
        "SUPERLU_DIST_OPENMP",
        "SUPERLUDIST_OpenMP",
        "OFF",
        "BOOL",
        "SuperLU_DIST OpenMP enabled",
        dependson="--superlu-dist",
    )

    # Magma
    group = parser.add_argument_group("MAGMA Options")

    add_arg(
        group, "--magma", "SUNDIALS_MAGMA", "ENABLE_MAGMA", "OFF", "BOOL", "SUNDIALS MAGMA support"
    )

    add_arg(
        group,
        "--magma-dir",
        "MAGMA_ROOT",
        "MAGMA_DIR",
        None,
        "PATH",
        "MAGMA install directory",
        dependson="--magma",
    )

    add_arg(
        group,
        "--magma-backends",
        "MAGAMA_BACKENDS",
        "SUNDIALS_MAGMA_BACKENDS",
        None,
        "STRING",
        "MAGMA backends",
        choices=["CUDA", "HIP"],
        dependson="--magma",
    )

    # ----------------
    # Other libraries
    # ----------------

    # hypre
    group = parser.add_argument_group("hypre Options")

    add_arg(
        group, "--hypre", "SUNDIALS_HYPRE", "ENABLE_HYPRE", "OFF", "BOOL", "SUNDIALS hypre support"
    )

    add_arg(
        group,
        "--hypre-incdir",
        "HYPRE_INCLUDE_DIR",
        "HYPRE_INCLUDE_DIR",
        None,
        "PATH",
        "Hypre include directory",
        dependson="--hypre",
    )

    add_arg(
        group,
        "--hypre-libdir",
        "HYPRE_LIBRARY_DIR",
        "HYPRE_LIBRARY_DIR",
        None,
        "PATH",
        "Hypre library directory",
        dependson="--hypre",
    )

    # PETSc
    group = parser.add_argument_group("PTESc Options")

    add_arg(
        group, "--petsc", "SUNDIALS_PETSC", "ENABLE_PETSC", "OFF", "BOOL", "SUNDIALS PETSc support"
    )

    add_arg(
        group,
        "--petsc-dir",
        "PETSC_ROOT",
        "PETSC_DIR",
        None,
        "PATH",
        "PETSc install directory",
        dependson="--petsc",
    )

    # Trilinos
    group = parser.add_argument_group("Trilinos Options")

    add_arg(
        group,
        "--trilinos",
        "SUNDIALS_TRILINOS",
        "ENABLE_TRILINOS",
        "OFF",
        "BOOL",
        "SUNDIALS Trilinos support",
    )

    add_arg(
        group,
        "--trilinos-dir",
        "TRILINOS_ROOT",
        "Trilinos_DIR",
        None,
        "PATH",
        "Trilinos install directory",
        dependson="--trilinos",
    )

    # XBraid
    group = parser.add_argument_group("XBraid Options")

    add_arg(
        group,
        "--xbraid",
        "SUNDIALS_XBRAID",
        "ENABLE_XBRAID",
        "OFF",
        "BOOL",
        "SUNDIALS XBraid support",
    )

    add_arg(
        group,
        "--xbraid-dir",
        "XBRAID_ROOT",
        "XBRAID_DIR",
        None,
        "PATH",
        "XBraid install directory",
        dependson="--xbraid",
    )

    # --------
    # Testing
    # --------

    group = parser.add_argument_group("Testing Options")

    # development tests
    add_arg(
        group,
        "--dev-tests",
        "SUNDIALS_TEST_ENABLE_DEV_TESTS",
        "SUNDIALS_TEST_ENABLE_DEV_TESTS",
        "OFF",
        "BOOL",
        "SUNDIALS development tests",
    )

    # unit tests
    add_arg(
        group,
        "--unit-tests",
        "SUNDIALS_TEST_ENABLE_UNIT_TESTS",
        "SUNDIALS_TEST_ENABLE_UNIT_TESTS",
        "OFF",
        "BOOL",
        "SUNDIALS unit tests",
    )

    add_arg(
        group,
        "--no-gtest",
        "SUNDIALS_TEST_ENABLE_GTEST",
        "SUNDIALS_TEST_ENABLE_GTEST",
        "ON",
        "BOOL",
        "SUNDIALS GTest unit tests",
    )

    # test output directory
    add_arg(
        group,
        "--test-output-dir",
        "SUNDIALS_TEST_OUTPUT_DIR",
        "SUNDIALS_TEST_OUTPUT_DIR",
        None,
        "PATH",
        "SUNDIALS test output directory",
    )

    # test answer directory
    add_arg(
        group,
        "--test-answer-dir",
        "SUNDIALS_TEST_ANSWER_DIR",
        "SUNDIALS_TEST_ANSWER_DIR",
        None,
        "PATH",
        "SUNDIALS test answer directory",
    )

    # test float comparison precision
    add_arg(
        group,
        "--test-float-precision",
        "SUNDIALS_TEST_FLOAT_PRECISION",
        "SUNDIALS_TEST_FLOAT_PRECISION",
        None,
        "STRING",
        "SUNDIALS test float comparison precision",
    )

    # test integer comparison precision
    add_arg(
        group,
        "--test-integer-precision",
        "SUNDIALS_TEST_INTEGER_PRECISION",
        "SUNDIALS_TEST_INTEGER_PRECISION",
        None,
        "STRING",
        "SUNDIALS test integer comparison precision",
    )

    add_arg(
        group,
        "--make-verbose",
        "CMAKE_VERBOSE_MAKEFILE",
        "CMAKE_VERBOSE_MAKEFILE",
        "OFF",
        "BOOL",
        "verbose make output",
    )

    # ---------------------
    # Parse and check args
    # ---------------------

    args = parser.parse_args()

    # output arg values from command line
    if args.debugscript:
        print_args(args)

    # read environment variables
    if args.readenv:
        read_env(args)

        # output arg values from env
        if args.debugscript:
            print_args(args)

    # ------------------
    # Create CMake file
    # ------------------

    with open(args.filename, "w") as fn:
        write_cmake(fn, args)


# -----------------------------------------------------------------------------
# Functions for reading environment variables
# -----------------------------------------------------------------------------


def read_env(args):
    """Set SUNDIALS options based on environment variables"""

    import os

    # get dictionary of input args
    args_dict = args.__dict__

    for a in args_dict:

        # skip non-cmake input args
        if type(args_dict[a]) is not dict:
            continue

        # don't overwrite options already set at command line
        value = args_dict[a]["value"]
        default = args_dict[a]["default"]

        if value != default:
            continue

        # check for environment variable and set value
        env_var = args_dict[a]["env_var"]

        if env_var is None:
            continue

        if env_var in os.environ:
            args_dict[a]["value"] = os.environ[env_var]


# -----------------------------------------------------------------------------
# Functions for adding script options
# -----------------------------------------------------------------------------


def add_arg(
    parser, arg, env_var, cmake_var, cmake_default, cmake_type, msg, choices=None, dependson=None
):
    """Add a command SUNDIALS option command line arg"""

    # Use underscores in the arg variable name
    arg_dest = arg[2:].replace("-", "_")
    help_msg = msg

    # Define function to create an argparse SUNDIALS option type
    arg_type = cmake_arg(
        env_var, cmake_var, cmake_default, cmake_type, msg, choices=choices, dependson=dependson
    )

    # Replace 'None' with a default string to ensure a dictionary is created
    # even when a command line input is not provided. This is ensures the
    # dictionary exists when reading variables from the environment.
    if cmake_default is None:
        cmake_default = "__default_none__"

    # Create command line arg
    parser.add_argument(arg, dest=arg_dest, type=arg_type, default=cmake_default, help=help_msg)


def cmake_arg(env_var, cmake_var, cmake_default, cmake_type, msg, choices=None, dependson=None):
    """Function factory for argparse SUNDIALS option type"""

    def _cmake_arg(str_var):
        """Create dictionary for SUNDIALS option"""

        import argparse

        # check if using None for the default value
        if str_var == "__default_none__":
            str_var = None

        # check for valid input options
        if cmake_type == "BOOL" and str_var not in ["ON", "OFF", None]:
            err_msg = "Invalid option value " + str_var + ". "
            err_msg += "Input value must be ON or OFF."
            raise argparse.ArgumentTypeError("Invalid Value for BOOL")

        if choices is not None and str_var is not None:
            raise_error = False
            if ";" in str_var:
                for s in str_var.split(";"):
                    if s not in choices:
                        raise_error = True
            else:
                if str_var not in choices:
                    raise_error = True

            if raise_error:
                err_msg = "Invalid option value " + str_var + ". "
                err_msg += "Input value must be "
                if len(choices) < 3:
                    err_msg += " or ".join(choices) + "."
                else:
                    err_msg += ", ".join(choices[:-1])
                    err_msg += ", or " + choices[-1] + "."
                raise argparse.ArgumentTypeError(err_msg)

        # create dictionary for SUNDIALS option
        cmake_dict = {}
        cmake_dict["env_var"] = env_var
        cmake_dict["cmake_var"] = cmake_var
        cmake_dict["default"] = cmake_default
        cmake_dict["cmake_type"] = cmake_type
        cmake_dict["msg"] = msg
        cmake_dict["value"] = str_var
        cmake_dict["depends_on"] = dependson

        return cmake_dict

    return _cmake_arg


# -----------------------------------------------------------------------------
# Functions for writing the CMake file
# -----------------------------------------------------------------------------


def write_cmake(fn, args):
    """Append commands to the CMake file"""

    # write output file comment header and (if necessary) script commands
    setup_file(fn, args.filename, args.filetype)

    # get dictionary of input args
    args_dict = args.__dict__

    for a in args_dict:

        # skip non-cmake input args
        if type(args_dict[a]) is not dict:
            continue

        # print(a, args_dict[a])

        # don't write output lines if using the default value
        value = args_dict[a]["value"]
        default = args_dict[a]["default"]

        if value is None or value == default:
            continue

        # don't write output if TPL is not enabled
        depends_on = args_dict[a]["depends_on"]

        if depends_on is not None:

            depends_on = depends_on[2:].replace("-", "_")
            depends_on_val = args_dict[depends_on]["value"]

            # print(depends_on, depends_on_val)

            if depends_on_val != "ON":
                continue

        # write CMake output
        cmake_var = args_dict[a]["cmake_var"]
        cmake_type = args_dict[a]["cmake_type"]
        cmake_msg = args_dict[a]["msg"]

        if args.filetype == "cache":
            cmd = f'set({cmake_var} "{value}" CACHE {cmake_type} ' f'"{cmake_msg}")\n'
        else:
            cmd = f' \\\n      -D {cmake_var}="{value}"'
        fn.write(cmd)


def setup_file(cmakefile, filename, filetype):
    """Setup the CMake file"""

    import os
    import stat

    if filetype == "cache":
        msg = (
            f"# CMake cache file for configuring SUNDIALS\n"
            f"#\n"
            f"# Move this file to your build directory and configure "
            f"SUNDIALS with the\n"
            f"# following command:\n"
            f"#   cmake <path to SUNDIALS source> -C {filename}\n"
        )
        cmakefile.write(msg)

        # update permissions to make sure the file is not executable
        NO_USER_EXE = ~stat.S_IXUSR
        NO_GROUP_EXE = ~stat.S_IXGRP
        NO_OTHER_EXE = ~stat.S_IXOTH
        NO_EXE = NO_USER_EXE & NO_GROUP_EXE & NO_OTHER_EXE

        st = os.stat(filename)
        os.chmod(filename, st.st_mode & NO_EXE)
    else:
        msg = (
            f"#!/bin/bash\n"
            f"# Script for configuring SUNDIALS\n"
            f"#\n"
            f"# Move this file to your build directory and configure "
            f"SUNDIALS with the\n"
            f"# following command:\n"
            f"#   ./{filename} <path to SUNDIALS source>\n"
            f'if [ "$#" -lt 1 ]; then\n'
            f'    echo "ERROR: Path to SUNDIALS source required"\n'
            f"    exit 1\n"
            f"fi\n"
            f"cmake $1"
        )
        cmakefile.write(msg)

        # update permissions to make sure the user can execute the script
        st = os.stat(filename)
        os.chmod(filename, st.st_mode | stat.S_IXUSR)


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------


def print_args(args):
    """Output command line input args"""

    for a in args.__dict__:
        if isinstance(args.__dict__[a], dict):
            print(a + ":")
            for key, value in args.__dict__[a].items():
                print("  ", key, "=", value)
        else:
            print(a, "=", args.__dict__[a])


# -----------------------------------------------------------------------------
# Run the main routine
# -----------------------------------------------------------------------------


if __name__ == "__main__":
    import sys

    sys.exit(main())
