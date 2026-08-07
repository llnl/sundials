#!/bin/bash
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
# SUNDIALS test environment setup script
#
# Usage: ./setup_env.sh <real type> <index size> <library type> <TPL status>
#                       <test type> [EXTRA ARGS]
#
# real type  = SUNDIALS real type to build/test with:
#                single   : single (32-bit) precision
#                double   : double (64-bit) precision
#                extended : extended (80-bit) precision
# index size = SUNDIALS index size to build/test with:
#                32       : 32-bit indices
#                64       : 64-bit indices
# lib type   = Which library type to test:
#                static   : only build static libraries
#                shared   : only build shared libraries
#                both     : build static and shared simultaneously
# TPL status = Enable/disable third party libraries:
#                ON       : All possible TPLs enabled
#                OFF      : No TPLs enabled
# test type  = Test type to run:
#                STD      : standard tests
#                DEV      : development tests
# EXTRA ARGS = Additional inputs passed to the environment script
# ------------------------------------------------------------------------------

echo "./setup_env.sh $*" | tee -a setup_env.log

# Check number of inputs
if [ "$#" -lt 5 ]; then
    echo "ERROR: missing required inputs"
    echo "  1) real type    : [single|double|extended]"
    echo "  2) index size   : [32|64]"
    echo "  3) library type : [static|shared|both]"
    echo "  4) TPL status   : [ON|OFF]"
    echo "  5) test type    : [STD|DEV]"
    return 1
fi

# These values should not be overridden by the environment script and may be
# referenced as appropriate to load TPLs
export SUNDIALS_PRECISION="$1"
export SUNDIALS_INDEX_SIZE="$2"
export SUNDIALS_LIBRARY_TYPE="$3"
export SUNDIALS_TPLS="$4"
export SUNDIALS_TEST_TYPE="$5"

# Remove parsed inputs, the remainder are passed to the environment script
shift 5

# Enable all SUNDIALS examples (the environment script may disable, CMake will
# disable if necessary depending on the configuration)
export SUNDIALS_EXAMPLES_C=ON
export SUNDIALS_EXAMPLES_CXX=ON
export SUNDIALS_EXAMPLES_F03=ON
export SUNDIALS_EXAMPLES_CUDA=ON

# Set the library type options (the environment script should not override)

case "${SUNDIALS_LIBRARY_TYPE}" in
    static)
        export SUNDIALS_STATIC_LIBRARIES=ON
        export SUNDIALS_SHARED_LIBRARIES=OFF
        ;;
    shared)
        export SUNDIALS_SHARED_LIBRARIES=ON
        export SUNDIALS_STATIC_LIBRARIES=OFF
        ;;
    both)
        export SUNDIALS_STATIC_LIBRARIES=ON
        export SUNDIALS_SHARED_LIBRARIES=ON
        ;;
    *)
        echo "ERROR: Unknown library type option: $SUNDIALS_LIBRARY_TYPE"
        return 1
        ;;
esac

# Set the test type options (the environment script should not override)

case "$SUNDIALS_TEST_TYPE" in
    DEV)
        export SUNDIALS_TEST_ENABLE_DEV_TESTS=ON
        export SUNDIALS_TEST_ENABLE_UNIT_TESTS=ON
        ;;
    STD|NONE)
        export SUNDIALS_TEST_ENABLE_DEV_TESTS=OFF
        export SUNDIALS_TEST_ENABLE_UNIT_TESTS=OFF
        ;;
    *)
        echo "ERROR: Unknown test type option: $SUNDIALS_TEST_TYPE"
        return 1
        ;;
esac

# Disable GTest in CI until we determine SEGFAULT cause in test_sundials_errors.c
export SUNDIALS_TEST_ENABLE_GTEST=OFF

# Build Type
export CMAKE_BUILD_TYPE="Debug"

# Disable compiler extensions by default. The user's environment script may
# override this setting if necessary.
export CMAKE_C_EXTENSIONS="OFF"
export CMAKE_CXX_EXTENSIONS="OFF"

# Enable compiler warnings (the environment script may disable)
export SUNDIALS_ENABLE_ALL_WARNINGS=ON
export SUNDIALS_ENABLE_WARNINGS_AS_ERRORS=ON

# Enable address sanitizer (environment script may disable)
# TODO(DJG): Always enable sanitizer not just when TPLs are OFF
if [[ "${SUNDIALS_TPLS}" == "OFF" ]]; then
    export SUNDIALS_ENABLE_ADDRESS_SANITIZER=ON
else
    export SUNDIALS_ENABLE_ADDRESS_SANITIZER=OFF
fi

# Call the environment setup script

if [ -n "$SUNDIALS_ENV_FILE" ]; then

    echo "Setting up environment with $SUNDIALS_ENV_FILE"
    # shellcheck source=/dev/null
    if ! source "$SUNDIALS_ENV_FILE" "$@"; then
        echo "ERROR: $SUNDIALS_ENV_FILE $* failed"
        return 1;
    fi

elif [ -f env/env.sh ]; then

    echo "Setting up environment with sundials/test/env/env.sh"
    # shellcheck source=/dev/null
    if ! source env/env.sh "$@"; then
        echo "ERROR: env/env.sh $* failed"
        return 1;
    fi

else

    echo "WARNING: No environment setup script found"

fi

# ensure TPLs are off when necessary

if [[ "${SUNDIALS_TPLS}" == "OFF" ]]; then

    # turn off profiling so we have a case that tests it
    export SUNDIALS_PROFILING=OFF

    # turn off logging so we have a case that tests it
    export SUNDIALS_LOGGING_LEVEL=0

    # threading
    export SUNDIALS_PTHREAD=OFF
    export SUNDIALS_OPENMP=OFF

    # mpi
    export SUNDIALS_MPI=OFF

    # gpu
    export SUNDIALS_CUDA=OFF
    export SUNDIALS_HIP=OFF
    export SUNDIALS_OPENMP_OFFLOAD=OFF

    # portability
    export SUNDIALS_KOKKOS=OFF
    export SUNDIALS_RAJA=OFF
    export SUNDIALS_SYCL=OFF

    # linear solvers
    export SUNDIALS_GINKGO=OFF
    export SUNDIALS_LAPACK=OFF
    export SUNDIALS_KLU=OFF
    export SUNDIALS_KOKKOS_KERNELS=OFF
    export SUNDIALS_SUPERLU_MT=OFF
    export SUNDIALS_SUPERLU_DIST=OFF
    export SUNDIALS_MAGMA=OFF

    # other libraries
    export SUNDIALS_HYPRE=OFF
    export SUNDIALS_PETSC=OFF
    export SUNDIALS_TRILINOS=OFF
    export SUNDIALS_XBRAID=OFF

    # fused kernels currently require CUDA or HIP
    export SUNDIALS_FUSED_KERNELS=OFF

fi

# Print relevant environment variables to the log file
for env_var in "${!CMAKE_@}"; do
    printf '%s=%s\n' "$env_var" "${!env_var}" >> setup_env.log
done
for env_var in "${!MPI@}"; do
    printf '%s=%s\n' "$env_var" "${!env_var}" >> setup_env.log
done
for env_var in "${!SUNDIALS_@}"; do
    printf '%s=%s\n' "$env_var" "${!env_var}" >> setup_env.log
done

# Check that only one of SUNDIALS_TEST_OUTPUT_DIR and SUNDIALS_TEST_ANSWER_DIR
# are set to ensure tests do not pass erroneously

if [ -n "${SUNDIALS_TEST_OUTPUT_DIR}" ] && [ -n "${SUNDIALS_TEST_ANSWER_DIR}" ]
then
    echo "ERROR: SUNDIALS_TEST_OUTPUT_DIR and SUNDIALS_TEST_ANSWER_DIR are set"
    echo "SUNDIALS_TEST_OUTPUT_DIR = ${SUNDIALS_TEST_OUTPUT_DIR}"
    echo "SUNDIALS_TEST_ANSWER_DIR = ${SUNDIALS_TEST_ANSWER_DIR}"
    return 1
fi

return 0
