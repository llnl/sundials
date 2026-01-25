#!/bin/bash
# ------------------------------------------------------------------------------
# Programmer(s): Cody J. Balos and David J. Gardner @ LLNL
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
# Script that sets up the default SUNDIALS testing environment.
#
# Usage: source docker.sh [compiler spec] [build type]
#
# Optional Inputs:
#   <compiler spec> = compiler to build SUNDIALS with
#   <build type>    = SUNDIALS build type:
#                       dbg : debug build
#                       opt : optimized build
# ------------------------------------------------------------------------------

echo "./docker.sh $*" | tee -a setup_env.log

# set defaults for optional inputs
compiler="gcc" # compiler spec
bldtype="dbg"  # build type dbg = debug or opt = optimized

if [ "$#" -gt 1 ]; then
    bldtype=$2
fi

# get compiler name and version from spec
compilername="${compiler%%@*}"
compilerversion="${compiler##*@}"

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

case "$compiler" in
    gcc);;
    *)
        echo "ERROR: Unknown compiler spec: $compiler"
        return 1
        ;;
esac

case "$bldtype" in
    dbg|opt) ;;
    *)
        echo "ERROR: Unknown build type: $bldtype"
        return 1
        ;;
esac

# ------------------------------------------------------------------------------
# Setup environment
# ------------------------------------------------------------------------------

# set file permissions (rwxrwxr-x)
umask 002

# path to shared installs
APPROOT=/opt

# setup the python environment
source ${APPROOT}/python-venv/sundocs/bin/activate

# ------------------------------------------------------------------------------
# Compilers and flags
# ------------------------------------------------------------------------------

if [ "$compilername" == "gcc" ]; then

    export CC=$(which gcc)
    export CXX=$(which g++)
    export FC=$(which gfortran)

    # optimization flags
    if [ "$bldtype" == "dbg" ]; then
        export CFLAGS="-g -O0"
        export CXXFLAGS="-g -O0"
        export FFLAGS="-g -O0 -fbounds-check"
        export CUDAFLAGS="-g -O0"
    else
        export CFLAGS="-g -O3"
        export CXXFLAGS="-g -O3"
        export FFLAGS="-g -O3"
        export CUDAFLAGS="-g -O3"
    fi

    # additional Fortran flags not currently added by ENABLE_ALL_WARNINGS
    export FFLAGS="${FFLAGS} -fcheck=all,no-pointer,no-recursion"
fi

# ------------------------------------------------------------------------------
# SUNDIALS Options
# ------------------------------------------------------------------------------

# Verbose build
export CMAKE_VERBOSE_MAKEFILE=OFF

# Number of build and test jobs
export SUNDIALS_BUILD_JOBS=4
export SUNDIALS_TEST_JOBS=1

# Sundials packages
export SUNDIALS_ARKODE=ON
export SUNDIALS_CVODE=ON
export SUNDIALS_CVODES=ON
export SUNDIALS_IDA=ON
export SUNDIALS_IDAS=ON
export SUNDIALS_KINSOL=ON

# Fortran interface status
if [ "$compilername" == "gcc" ]; then
    if [[ ("$SUNDIALS_PRECISION" == "double") ]]; then
        export SUNDIALS_FMOD_INTERFACE=ON
    else
        export SUNDIALS_FMOD_INTERFACE=OFF
    fi
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

# Answer files
if [ -z "${SUNDIALS_TEST_ANSWER_DIR}" ]; then
    export SUNDIALS_TEST_ANSWER_DIR="${PWD}/answers/linux-ubuntu20.04-x86_64/gcc-9.4.0/${SUNDIALS_PRECISION}"
fi

# Uncomment to override the default output file comparison precisions. The float
# precision is number of digits to compare (0 = all digits) and the integer
# precision is allowed percentage difference (0 = no difference).
if [ "$SUNDIALS_PRECISION" == "extended" ]; then
    export SUNDIALS_TEST_FLOAT_PRECISION=7
    export SUNDIALS_TEST_INTEGER_PRECISION=10
elif [ "$SUNDIALS_PRECISION" == "double" ]; then
    export SUNDIALS_TEST_FLOAT_PRECISION=5
    export SUNDIALS_TEST_INTEGER_PRECISION=5
else # single
    export SUNDIALS_TEST_FLOAT_PRECISION=3
    export SUNDIALS_TEST_INTEGER_PRECISION=10
fi

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

# ---
# MPI
# ---

MPI_ROOT=/opt/view

export SUNDIALS_MPI=ON
export MPICC="${MPI_ROOT}/bin/mpicc"
export MPICXX="${MPI_ROOT}/bin/mpicxx"
export MPIFC="${MPI_ROOT}/bin/mpifort"
export MPIEXEC="${MPI_ROOT}/bin/mpirun"
export MPIEXEC_PREFLAGS="--oversubscribe"

# -------------
# LAPACK / BLAS
# -------------

if [ "$SUNDIALS_PRECISION" == "single" ]; then
    export SUNDIALS_LAPACK=ON
    export LAPACK_ROOT=/opt/view
    export LAPACK_LIBRARIES="${LAPACK_ROOT}/lib/libopenblas.so"
elif [ "$SUNDIALS_PRECISION" == "double" ]; then
    export SUNDIALS_LAPACK=ON
    export LAPACK_ROOT=/opt/view
    export LAPACK_LIBRARIES="${LAPACK_ROOT}/lib/libopenblas.so"
else
    export SUNDIALS_LAPACK=OFF
    unset LAPACK_LIBRARIES
fi

# ---
# KLU
# ---

if [ "$SUNDIALS_PRECISION" == "double" ]; then
    export SUNDIALS_KLU=ON
    export SUITE_SPARSE_ROOT=/opt/view
    export SUITE_SPARSE_INCLUDE_DIR="${SUITE_SPARSE_ROOT}/include"
    export SUITE_SPARSE_LIBRARY_DIR="${SUITE_SPARSE_ROOT}/lib"
else
    export SUNDIALS_KLU=OFF
    unset SUITE_SPARSE_ROOT
fi

# ----------
# SuperLU_MT
# ----------

if [ "$SUNDIALS_PRECISION" == "single" ]; then
    export SUNDIALS_SUPERLU_MT=ON
    export SUPERLU_MT_ROOT=/opt/view/superlu-mt
    export SUPERLU_MT_INCLUDE_DIR="${SUPERLU_MT_ROOT}/include"
    export SUPERLU_MT_LIBRARY_DIR="${SUPERLU_MT_ROOT}/lib"
    export SUPERLU_MT_LIBRARIES="${SUPERLU_MT_ROOT}/lib/libblas_PTHREAD.a;${SUPERLU_MT_ROOT}/lib/libsuperlu_mt_PTHREAD.a"
    export SUPERLU_MT_THREAD_TYPE="PTHREAD"
elif [ "$SUNDIALS_PRECISION" == "double" ]; then
    export SUNDIALS_SUPERLU_MT=ON
    export SUPERLU_MT_ROOT=/opt/view/superlu-mt
    export SUPERLU_MT_INCLUDE_DIR="${SUPERLU_MT_ROOT}/include"
    export SUPERLU_MT_LIBRARY_DIR="${SUPERLU_MT_ROOT}/lib"
    export SUPERLU_MT_LIBRARIES="${SUPERLU_MT_ROOT}/lib/libblas_PTHREAD.a;${SUPERLU_MT_ROOT}/lib/libsuperlu_mt_PTHREAD.a"
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
    export SUPERLU_DIST_ROOT=/opt/view/superlu-dist
    export SUPERLU_DIST_INCLUDE_DIR="${SUPERLU_DIST_ROOT}/include"
    export SUPERLU_DIST_LIBRARY_DIR="${SUPERLU_DIST_ROOT}/lib"

    # build with netlib blas/lapack if using 64-bit indices so that we can
    # use 64-bit openblas for other TPLs
    export BLAS_ROOT=/opt/view
    if [ "$SUNDIALS_INDEX_SIZE" == "64" ]; then
        if [ -f "${BLAS_ROOT}/lib/libblas.so" ]; then
            export BLAS_LIBRARIES="${BLAS_ROOT}/lib/libblas.so;${BLAS_ROOT}/lib/liblapack.so"
        fi
    else
        export BLAS_LIBRARIES="${BLAS_ROOT}/lib/libopenblas.so"
    fi

    # PARMETIS
    export PARMETIS_ROOT=/opt/view
    export PARMETIS_LIBRARIES="${PARMETIS_ROOT}/lib/libparmetis.so"

    # METIS
    export METIS_ROOT=/opt/view
    export METIS_LIBRARIES="${METIS_ROOT}/lib/libmetis.so"

    export SUPERLU_DIST_LIBRARIES="${BLAS_LIBRARIES};${PARMETIS_LIBRARIES};${METIS_LIBRARIES};${SUPERLU_DIST_ROOT}/lib/libsuperlu_dist.a"
    export SUPERLU_DIST_OPENMP=OFF

    # if BLAS wasn't found, then dont build SuperLU_DIST
    if [ -z "$BLAS_LIBRARIES" ]; then
        export SUNDIALS_SUPERLU_DIST=OFF
    else
        export SUNDIALS_SUPERLU_DIST=ON
    fi
else
    export SUNDIALS_SUPERLU_DIST=OFF
    unset SUPERLU_DIST_INCLUDE_DIR
    unset SUPERLU_DIST_LIBRARY_DIR
    unset SUPERLU_DIST_LIBRARIES
    unset SUPERLU_DIST_OPENMP
fi

# -----
# hypre
# -----

if [ "$SUNDIALS_PRECISION" == "double" ]; then
    export SUNDIALS_HYPRE=ON
    export HYPRE_ROOT=/opt/view
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
    export PETSC_ROOT=/opt/view
else
    export SUNDIALS_PETSC=OFF
    unset PETSC_ROOT
fi

# --------
# trilinos
# --------

if [ "$SUNDIALS_PRECISION" == "double" ] && [ "$SUNDIALS_INDEX_SIZE" == "32" ]; then
    export SUNDIALS_TRILINOS=OFF
    export TRILINOS_ROOT=/opt/view
else
    export SUNDIALS_TRILINOS=OFF
    unset TRILINOS_ROOT
fi

# ------
# xbraid
# ------

if [ "$SUNDIALS_PRECISION" == "double" ] && [ "$SUNDIALS_INDEX_SIZE" == "32" ]; then
    export SUNDIALS_XBRAID=ON
    export XBRAID_ROOT=/opt/view
else
    export SUNDIALS_XBRAID=OFF
    unset XBRAID_ROOT
fi
