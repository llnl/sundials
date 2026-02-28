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
# Remove carriage returns (^M) from all files in a given directory recursively
# ------------------------------------------------------------------------------

# check number of inputs
if [ "$#" -lt 1 ]; then
    echo "ERROR: Search directory path required"
    exit 1
fi

# wrapper for that appends the correct flags for editing inplace with OS X or
# Linux sed implementations
sedi() {
    case $(uname) in
        Darwin*) sedi=('-i' '') ;;
        *) sedi='-i' ;;
    esac
    sed "${sedi[@]}" "$@"
}

# remove carriage returns from all files in the given directory recursively
for f in $(grep -Ilr $'\r' $1); do
    echo "Updating: $f"
    sedi $'s/\r//g' $f
done

echo "Done"
