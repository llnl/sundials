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
# code, cmake-format for CMake files, and ruff for Python code.
#
# Usage:
#    ./format.sh [--clang-format|--fprettify|--cmake|--ruff] \
#               [paths to directories or files to format]
#
# We require clang-format 17.0.4. Other versions may produce different styles!
# ---------------------------------------------------------------------------------

print_usage() {
    cat <<EOF
Usage: format.sh [formatter] [paths]

Formatters:
  --clang-format  Format C/C++/CUDA files with clang-format.
  --fprettify     Format Fortran files with fprettify.
  --cmake         Format CMake files with cmake-format.
  --ruff           Format Python files with ruff.

If no formatter is specified, all formatters are run. If a formatter is
specified without any paths, the current directory is formatted.
EOF
}

formatter=all

while [[ $# -gt 0 ]]; do
    case "$1" in
        --help|-h)
            print_usage
            exit 0
            ;;
        --clang-format|--fprettify|--cmake|--cmake-format|--ruff)
            if [[ "$formatter" != all ]]; then
                echo "ERROR: Specify at most one formatter" >&2
                exit 1
            fi
            formatter="${1#--}"
            shift
            ;;
        --)
            shift
            break
            ;;
        -* )
            echo "ERROR: Unknown option: $1" >&2
            print_usage >&2
            exit 1
            ;;
        *)
            break
            ;;
    esac
done

if [ $# -lt 1 ]; then
    if [[ "$formatter" == all ]]; then
        echo "ERROR: At least one path to format required" >&2
        exit 1
    fi
    paths=(.)
else
    paths=( "$@" )
fi

if [[ "$formatter" == all || "$formatter" == clang-format ]]; then
    find "${paths[@]}" \( -iname '*.h' -o -iname '*.hpp' -o \
         -iname '*.c' -o -iname '*.cpp' -o \
         -iname '*.cuh' -o -iname '*.cu' \) -print | grep -v fmod | \
         xargs -r clang-format -i
fi

if [[ "$formatter" == all || "$formatter" == fprettify ]]; then
    find "${paths[@]}" -iname '*.f90' -print | grep -v fmod | \
         xargs -r fprettify --indent 2 --enable-replacements --c-relations
fi

if [[ "$formatter" == all || "$formatter" == cmake || "$formatter" == cmake-format ]]; then
    find "${paths[@]}" \( -iname '*.cmake' -o -iname 'CMakeLists.txt' \) \
         -exec cmake-format -i {} ';'
fi

if [[ "$formatter" == all || "$formatter" == ruff ]]; then
    ruff format "${paths[@]}"
fi
