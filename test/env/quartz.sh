#!/bin/bash
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
# Script that sets up the environment for SUNDIALS testing on Quartz.
#
# Usage: source env.quartz.sh <compiler spec> <build type>
#
# Optional Inputs:
#   <compiler spec> = Compiler to build sundials with:
#                       e.g., gcc@8.1.0, intel@19.0.4, clang@9.0.0, etc.
#   <build type>    = SUNDIALS build type:
#                       dbg : debug build
#                       opt : optimized build
# ------------------------------------------------------------------------------

echo "./quartz.sh $*" | tee -a setup_env.log

# set defaults for optional inputs
compiler="gcc@8.1.0" # compiler spec
bldtype="dbg"        # build type dbg = debug or opt = optimized

# set optional inputs if provided
if [ "$#" -gt 0 ]; then
    compiler=$1
fi

if [ "$#" -gt 1 ]; then
    bldtype=$2
fi

# ------------------------------------------------------------------------------
# Setup environment
# ------------------------------------------------------------------------------

# get compiler name and version from spec
compilername="${compiler%%@*}"
compilerversion="${compiler##*@}"

# load default cmake module
module load cmake

# load compiler module
case "$compilername" in
    clang)
        module load clang/${compilerversion}
        if [ $? -ne 0 ]; then return 1; fi
        export CC=$(which clang)
        export CXX=$(which clang++)
        export FC=""
        FORTRAN_STATUS=OFF
        ;;
    gcc)
        module load gcc/${compilerversion}
        if [ $? -ne 0 ]; then return 1; fi
        export CC=$(which gcc)
        export CXX=$(which g++)
        export FC=$(which gfortran)
        FORTRAN_STATUS=ON
        ;;
    intel)
        module load intel/${compilerversion}
        if [ $? -ne 0 ]; then return 1; fi
        export CC=$(which icc)
        export CXX=$(which icpc)
        export FC=$(which ifort)
        FORTRAN_STATUS=ON
        ;;
    pgi)
        module load pgi/${compilerversion}
        if [ $? -ne 0 ]; then return 1; fi
        export CC=$(which pgcc)
        export CXX=$(which pgc++)
        export FC=$(which pgfortran)
        FORTRAN_STATUS=ON
        ;;
    *)
        echo "ERROR: Invalid compiler option: $compiler"
        return 1
        ;;
esac

# set compiler flags
if [ "$bldtype" == "dbg" ]; then
    export CFLAGS="-g -O0"
    export CXXFLAGS="-g -O0"
    export FFLAGS="-g -O0"
elif [ "$bldtype" == "opt" ]; then
    export CFLAGS="-g -O3"
    export CXXFLAGS="-g -O3"
    export FFLAGS="-g -O3"
else
    echo "ERROR: Invalid build type: $bldtype"
    return 1
fi

# Fortran settings
if [[ ("$SUNDIALS_PRECISION" == "double") ]]; then
    export SUNDIALS_FMOD_INTERFACE=${FORTRAN_STATUS}
else
    export SUNDIALS_FMOD_INTERFACE=OFF
fi

# set MPI compiler wrapper
export SUNDIALS_MPI=ON
export MPICC=$(which mpicc)
export MPICXX=$(which mpicxx)
export MPIFC=$(which mpifort)
export MPIEXEC=$(which srun)

# PThread settings
export SUNDIALS_PTHREAD=ON

# OpenMP settings
export SUNDIALS_OPENMP=ON
export OMP_NUM_THREADS=16
