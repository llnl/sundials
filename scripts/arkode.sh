#!/bin/bash
# ------------------------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds, David J. Gardner, Cody J. Balos @ LLNL
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
# Script to add ARKODE files to a SUNDIALS tar-file.
# ------------------------------------------------------------------------------

set -e
set -o pipefail

tarfile=$1
distrobase=$2
doc=$3

# all remaining inputs are for tar command
shift 3
tar=$*

echo "   --- Add arkode module to $tarfile"

if [ $doc = "T" ]; then
    $tar $tarfile $distrobase/doc/arkode/ark_guide.pdf
    $tar $tarfile $distrobase/doc/arkode/ark_examples.pdf
fi

echo "   --- Add arkode include files to $tarfile"
$tar $tarfile $distrobase/include/arkode

echo "   --- Add arkode source files to $tarfile"
$tar $tarfile $distrobase/src/arkode

echo "   --- Add arkode examples to $tarfile"
$tar $tarfile $distrobase/examples/arkode

echo "   --- Add arkode unit tests to $tarfile"
$tar $tarfile $distrobase/test/unit_tests/arkode
