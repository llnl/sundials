#!/bin/bash
# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos and David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------
# Script that sets up the default SUNDIALS testing environment.
# ------------------------------------------------------------------------------

echo "./jenkins.sh $*" | tee -a setup_env.log

# ------------------------------------------------------------------------------
# Check input values
# ------------------------------------------------------------------------------

case "$SUNDIALS_PRECISION" in
    single|double|extended) ;;
    *)
        echo "ERROR: Unknown real type option: $SUNDIALS_PRECISION"
        return 1
        ;;
esac

case "$SUNDIALS_INDEX_SIZE" in
    32|64) ;;
    *)
        echo "ERROR: Unknown index size option: $SUNDIALS_INDEX_SIZE"
        return 1
        ;;
esac

# ------------------------------------------------------------------------------
# Setup environment
# ------------------------------------------------------------------------------

# set file permissions (rwxrwxr-x)
umask 002

# setup the python environment
source /usr/local/suntest/pyevn/sundials/bin/activate

# setup spack
export SPACK_ROOT=/usr/local/suntest/spack
source ${SPACK_ROOT}/share/spack/setup-env.sh

# add CUDA
export CUDAARCHS=60

if [[ ":${PATH}:" != *":/usr/local/cuda/bin:"* ]]; then
    export PATH="/usr/local/cuda/bin${PATH:+:${PATH}}"
fi

if [[ ":${LD_LIBRARY_PATH}:" != *":/usr/local/cuda/lib64:"* ]]; then
    export LD_LIBRARY_PATH="/usr/local/cuda/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
fi

# load CMake
spack load cmake@3.18.6

# ------------------------------------------------------------------------------
# Compilers and flags
# ------------------------------------------------------------------------------

export CC="$(which gcc)"
export CXX="$(which g++)"
export FC="$(which gfortran)"

# disable optimization
export CFLAGS="-O0"
export CXXFLAGS="-O0"
export FFLAGS="-O0"
export CUDAFLAGS="-O0"

# ------------------------------------------------------------------------------
# SUNDIALS Options
# ------------------------------------------------------------------------------

# Verbose build
export CMAKE_VERBOSE_MAKEFILE=OFF

export CMAKE_BUILD_TYPE=Debug

# Number of build and test jobs
export SUNDIALS_BUILD_JOBS=6
export SUNDIALS_TEST_JOBS=1

# Sundials packages
export SUNDIALS_ARKODE=ON
export SUNDIALS_CVODE=ON
export SUNDIALS_CVODES=ON
export SUNDIALS_IDA=ON
export SUNDIALS_IDAS=ON
export SUNDIALS_KINSOL=ON

# Fortran interface status
if [[ ("$SUNDIALS_PRECISION" == "double") ]]; then
    export SUNDIALS_FMOD_INTERFACE=ON
else
    export SUNDIALS_FMOD_INTERFACE=OFF
fi

# Sundials benchmarks
export SUNDIALS_BENCHMARKS=ON

# Sundials monitoring
export SUNDIALS_MONITORING=ON

# Sundials profiling
export SUNDIALS_PROFILING=ON

# Sundials logging
export SUNDIALS_LOGGING_LEVEL=3

# Uncomment to override the default output file comparison precisions. The float
# precision is number of digits to compare (0 = all digits) and the integer
# precision is allowed percentage difference (0 = no difference).
export SUNDIALS_TEST_FLOAT_PRECISION=0
export SUNDIALS_TEST_INTEGER_PRECISION=0

# ------------------------------------------------------------------------------
# Third party libraries
# ------------------------------------------------------------------------------

# -------
# PThread
# -------

export SUNDIALS_PTHREAD=ON

# ------
# OpenMP
# ------

export SUNDIALS_OPENMP=ON
export OMP_NUM_THREADS=4
export OMP_PROC_BIND=false

# ---------------------
# OpenMP Device Offload
# ---------------------

export SUNDIALS_OPENMP_OFFLOAD=OFF

# ----
# CUDA
# ----

if [ "$SUNDIALS_PRECISION" != "extended" ]; then
    export SUNDIALS_CUDA=ON
    export SUNDIALS_FUSED_KERNELS=ON
else
    export SUNDIALS_CUDA=OFF
    export SUNDIALS_FUSED_KERNELS=OFF
fi

# ---
# MPI
# ---

MPI_ROOT="$(spack location -i openmpi@5.0.5)"

export SUNDIALS_MPI=ON
export MPICC="${MPI_ROOT}/bin/mpicc"
export MPICXX="${MPI_ROOT}/bin/mpicxx"
export MPIFC="${MPI_ROOT}/bin/mpifort"
export MPIEXEC="${MPI_ROOT}/bin/mpirun"

# -------------
# LAPACK / BLAS
# -------------

if [ "$SUNDIALS_PRECISION" != "extended" ]; then
    export SUNDIALS_LAPACK=ON
    if [ "$SUNDIALS_INDEX_SIZE" == "32" ]; then
        LAPACK_ROOT="$(spack location -i openblas@0.3.27 ~ilp64)"
    else
        LAPACK_ROOT="$(spack location -i openblas@0.3.27 +ilp64)"
    fi
    export LAPACK_LIBRARIES="${LAPACK_ROOT}/lib/libopenblas.so"
    export BLAS_LIBRARIES="${LAPACK_LIBRARIES}"
else
    export SUNDIALS_LAPACK=OFF
    unset LAPACK_LIBRARIES
    unset BLAS_LIBRARIES
fi

# ---
# KLU
# ---

if [ "$SUNDIALS_PRECISION" == "double" ]; then
    export SUNDIALS_KLU=ON
    if [ "$SUNDIALS_INDEX_SIZE" == "32" ]; then
        SUITE_SPARSE_ROOT="$(spack location -i suite-sparse@7.7.0 ^openblas ~ilp64)"
    else
        SUITE_SPARSE_ROOT="$(spack location -i suite-sparse@7.7.0 ^openblas +ilp64)"
    fi
    export SUITE_SPARSE_INCLUDE_DIR="${SUITE_SPARSE_ROOT}/include"
    export SUITE_SPARSE_LIBRARY_DIR="${SUITE_SPARSE_ROOT}/lib"
else
    export SUNDIALS_KLU=OFF
    unset SUITE_SPARSE_ROOT
fi

# ------
# Ginkgo
# ------

if [ "$SUNDIALS_PRECISION" != "extended" ]; then
    if [ "$SUNDIALS_CUDA" == "ON" ]; then
        export SUNDIALS_GINKGO=ON
        export GINKGO_ROOT="$(spack location -i ginkgo@master +cuda)"
        export GINKGO_BACKENDS="REF;OMP;CUDA"
    else
        export SUNDIALS_GINKGO=ON
        export GINKGO_ROOT="$(spack location -i ginkgo@master ~cuda)"
        export GINKGO_BACKENDS="REF;OMP"
    fi
else
    export SUNDIALS_GINKGO=OFF
    unset GINKGO_ROOT
    unset GINKGO_BACKENDS
fi

# ------
# Kokkos
# ------

if [ "$SUNDIALS_PRECISION" == "double" ]; then
    export SUNDIALS_KOKKOS=ON
    export KOKKOS_ROOT="$(spack location -i kokkos@4.3.01 ~cuda)"
else
    export SUNDIALS_KOKKOS=OFF
    unset KOKKOS_ROOT
fi

# --------------
# Kokkos-Kernels
# --------------

if [ "$SUNDIALS_PRECISION" == "double" ]; then
    export SUNDIALS_KOKKOS_KERNELS=ON
    export KOKKOS_KERNELS_ROOT="$(spack location -i kokkos-kernels@4.3.01 ~cuda)"
else
    export SUNDIALS_KOKKOS_KERNELS=OFF
    unset KOKKOS_KERNELS_ROOT
fi

# -----
# MAGMA
# -----

if [ "$SUNDIALS_PRECISION" != "extended" ] && \
    [ "$SUNDIALS_INDEX_SIZE" == "32" ] && \
    [ "$SUNDIALS_CUDA" == "ON" ]; then
    export SUNDIALS_MAGMA=ON
    export MAGMA_ROOT="$(spack location -i magma@2.8.0 +cuda)"
    export MAGMA_BACKENDS="CUDA"
else
    export SUNDIALS_MAGMA=OFF
    unset MAGMA_ROOT
    unset MAGMA_BACKENDS
fi

# ----------
# SuperLU_MT
# ----------

if [ "$SUNDIALS_PRECISION" != "extended" ]; then
    export SUNDIALS_SUPERLU_MT=ON
    # Using @master (sha 9e23fe72652afc28c97829e69e7c6966050541a7) as it
    # additional fixes necessary for building with newer versions of GCC
    if [ "$SUNDIALS_INDEX_SIZE" == "32" ]; then
        SUPERLU_MT_ROOT="$(spack location -i superlu-mt@master ~int64)"
    else
        SUPERLU_MT_ROOT="$(spack location -i superlu-mt@master +int64)"
    fi
    export SUPERLU_MT_INCLUDE_DIR="${SUPERLU_MT_ROOT}/include"
    export SUPERLU_MT_LIBRARY_DIR="${SUPERLU_MT_ROOT}/lib"
    export SUPERLU_MT_LIBRARIES="${SUPERLU_MT_ROOT}/lib/libblas_PTHREAD.a"
    export SUPERLU_MT_THREAD_TYPE="PTHREAD"
else
    export SUNDIALS_SUPERLU_MT=OFF
    unset SUPERLU_MT_ROOT
    unset SUPERLU_MT_INCLUDE_DIR
    unset SUPERLU_MT_LIBRARY_DIR
    unset SUPERLU_MT_LIBRARIES
    unset SUPERLU_MT_THREAD_TYPE
fi

# ------------
# SuperLU_DIST
# ------------

if [ "$SUNDIALS_PRECISION" == "double" ]; then
    export SUNDIALS_SUPERLU_DIST=ON
    if [ "$SUNDIALS_INDEX_SIZE" == "32" ]; then
        export SUPERLU_DIST_ROOT="$(spack location -i superlu-dist@8.2.1 ~int64 ~cuda)"
    else
        export SUPERLU_DIST_ROOT="$(spack location -i superlu-dist@8.2.1 +int64 ~cuda)"
    fi
    export SUPERLU_DIST_OPENMP=OFF
else
    export SUNDIALS_SUPERLU_DIST=OFF
    unset SUPERLU_DIST_ROOT
    unset SUPERLU_DIST_OPENMP
fi

# -----
# hypre
# -----

if [ "$SUNDIALS_PRECISION" == "double" ]; then
    export SUNDIALS_HYPRE=ON
    if [ "$SUNDIALS_INDEX_SIZE" == "32" ]; then
        HYPRE_ROOT="$(spack location -i hypre@2.31.0 ~int64 ~cuda)"
    else
        HYPRE_ROOT="$(spack location -i hypre@2.31.0 +int64 ~cuda)"
    fi
    export HYPRE_INCLUDE_DIR="${HYPRE_ROOT}/include"
    export HYPRE_LIBRARY_DIR="${HYPRE_ROOT}/lib"
else
    export SUNDIALS_HYPRE=OFF
    unset HYPRE_INCLUDE_DIR
    unset HYPRE_LIBRARY_DIR
fi

# -----
# petsc
# -----

if [ "$SUNDIALS_PRECISION" == "double" ]; then
    export SUNDIALS_PETSC=ON
    if [ "$SUNDIALS_INDEX_SIZE" == "32" ]; then
        PETSC_ROOT="$(spack location -i petsc@3.21.4 +double ~int64 ~cuda)"
    else
        PETSC_ROOT="$(spack location -i petsc@3.21.4 +double +int64 ~cuda)"
    fi
    export PETSC_ROOT
else
    export SUNDIALS_PETSC=OFF
    unset PETSC_ROOT
fi

# --------
# trilinos
# --------

if [ "$SUNDIALS_PRECISION" == "double" ] && [ "$SUNDIALS_INDEX_SIZE" == "32" ]; then
    export SUNDIALS_TRILINOS=ON
    if [ "$SUNDIALS_INDEX_SIZE" == "32" ]; then
        TRILINOS_ROOT="$(spack location -i trilinos@16.0.0 gotype=int ~cuda)"
    else
        TRILINOS_ROOT="$(spack location -i trilinos@16.0.0 gotype=long_long ~cuda)"
    fi
    export TRILINOS_ROOT
else
    export SUNDIALS_TRILINOS=OFF
    unset TRILINOS_ROOT
fi

# ----
# raja
# ----

if [ "$SUNDIALS_PRECISION" == "double" ]; then
    export SUNDIALS_RAJA=ON
    export RAJA_ROOT="$(spack location -i raja@2024.02.2 +cuda)"
    export RAJA_BACKENDS="CUDA"
    # RAJA does not find camp on its own?
    export camp_ROOT="$(spack location -i camp@2024.02.1 +cuda)"
else
    export SUNDIALS_RAJA=OFF
    unset RAJA_ROOT
    unset RAJA_BACKENDS
    unset camp_ROOT
fi

# ------
# xbraid
# ------

if [ "$SUNDIALS_PRECISION" == "double" ] && [ "$SUNDIALS_INDEX_SIZE" == "32" ]; then
    export SUNDIALS_XBRAID=ON
    export XBRAID_ROOT="$(spack location -i xbraid@3.0.0)"
else
    export SUNDIALS_XBRAID=OFF
    unset XBRAID_ROOT
fi
