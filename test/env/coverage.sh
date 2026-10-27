#!/bin/bash
# ------------------------------------------------------------------------------
# Programmer(s): SUNDIALS contributors
# ------------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# ------------------------------------------------------------------------------
# Coverage test environment based on the standard CI environment.

# shellcheck source=env/docker.sh
source "$(dirname "${BASH_SOURCE[0]}")/docker.sh" "$@" || return 1

# Let CMake add coverage instrumentation and provide the `codecov` target.
export SUNDIALS_ENABLE_CODE_COVERAGE=ON

# setup_env.sh sets Debug by default before sourcing this file. Allow the
# workflow (or a developer running test_driver.sh locally) to override it.
if [[ -n "${SUNDIALS_COVERAGE_BUILD_TYPE:-}" ]]; then
    case "${SUNDIALS_COVERAGE_BUILD_TYPE}" in
        Debug|RelWithDebInfo)
            export CMAKE_BUILD_TYPE="${SUNDIALS_COVERAGE_BUILD_TYPE}"
            ;;
        *)
            echo "ERROR: Unknown coverage build type ${SUNDIALS_COVERAGE_BUILD_TYPE}"
            return 1
            ;;
    esac
fi
