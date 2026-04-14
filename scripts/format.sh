#!/bin/bash
# ---------------------------------------------------------------------------------
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
# ---------------------------------------------------------------------------------
# This script will use clang-format to format C/C++ code, fprettify for Fortran
# code, cmake-format for CMake files, and black for Python code.
#
# Usage:
#    ./format.sh <paths to directories or files to format>
#
# We require clang-format 17.0.4. Other versions may produce different styles!
# ---------------------------------------------------------------------------------

if [ $# -lt 1 ]; then
    echo "ERROR: At least one path to format required"
    exit 1
fi

paths=( "$@" )

find "${paths[@]}" -iname '*.h' -o -iname '*.hpp' -o \
     -iname '*.c' -o -iname '*.cpp' -o \
     -iname '*.cuh' -o -iname '*.cu' | grep -v fmod | xargs clang-format -i

find "${paths[@]}" -iname '*.f90' | grep -v fmod | xargs fprettify --indent 2 --enable-replacements --c-relations

find "${paths[@]}" \( -iname '*.cmake' -o -iname 'CMakeLists.txt' \) \
     -exec cmake-format -i {} ';'

find "${paths[@]}" -iname '*.py' -exec black {} ';'
