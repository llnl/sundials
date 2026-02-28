# ---------------------------------------------------------------
# Programmer(s): Cody Balos @ LLNL
# ---------------------------------------------------------------
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
# ---------------------------------------------------------------
# SuperLUDIST find module that creates an imported target for
# SuperLU_DIST. The target is SUNDIALS::SUPERLUDIST.
#
# The module has two modes: one 'automatic' mode that uses
# pkg-conf to locate and determine the SuperLU DIST
# libraries and includes, and one 'manual' mode which
# requires the libraries and includes to be specified.
#
# If SUPERLUDIST_LIBRARIES is set, then the 'manual' mode is
# chosen and the SUPERLUDIST_INCLUDE_DIRS variable must also
# be set.
#
# The variable SUPERLUDIST_DIR can be used to control where
# the module looks for the library in 'automatic' mode.
#
# This module also defines additional variables:
#
#   SUPERLUDIST_FOUND      - the SuperLU_DIST libraries were found
#   SUPERLUDIST_INDEX_SIZE - the bit width of indices in SUPERLUDIST
#   SUPERLUDIST_CUDA       - SuperLU_DIST has CUDA support
#   SUPERLUDIST_ROCM       - SuperLU_DIST has ROCm support
# ---------------------------------------------------------------

if(NOT SUPERLUDIST_LINK_LIBRARIES AND SUPERLUDIST_LIBRARIES)
  set(SUPERLUDIST_LINK_LIBRARIES
      "${SUPERLUDIST_LIBRARIES}"
      CACHE INTERNAL "")
elseif(NOT SUPERLUDIST_LINK_LIBRARIES)
  find_package(PkgConfig REQUIRED)
  list(APPEND CMAKE_PREFIX_PATH "${SUPERLUDIST_DIR}")
  if(DEFINED SUPERLUDIST_FIND_VERSION)
    if(SUPERLUDIST_FIND_VERSION_EXACT)
      set(_pkg_version_spec "=${SUPERLUDIST_FIND_VERSION}")
    else()
      set(_pkg_version_spec ">=${SUPERLUDIST_FIND_VERSION}")
    endif()
  endif()
  pkg_search_module(SUPERLUDIST REQUIRED "superlu_dist${_pkg_version_spec}")
  set(SUPERLUDIST_LINK_LIBRARIES
      "${SUPERLUDIST_LINK_LIBRARIES}"
      CACHE INTERNAL "")
  set(SUPERLUDIST_INCLUDE_DIRS
      "${SUPERLUDIST_INCLUDE_DIRS}"
      CACHE INTERNAL "")
endif()

# find the library configuration file
set(SUPERLUDIST_CUDA
    FALSE
    CACHE BOOL "SuperLU DIST was built with CUDA support")
set(SUPERLUDIST_ROCM
    FALSE
    CACHE BOOL "SuperLU DIST was built with ROCm support")
if(SUPERLUDIST_INCLUDE_DIRS)
  find_file(SUPERLUDIST_CONFIG_PATH superlu_dist_config.h
            PATHS "${SUPERLUDIST_INCLUDE_DIRS}")
  mark_as_advanced(FORCE SUPERLUDIST_CONFIG_PATH)
  if(SUPERLUDIST_VERSION VERSION_GREATER_EQUAL "8.0.0")
    file(STRINGS "${SUPERLUDIST_CONFIG_PATH}" _index_size_64
         REGEX "#define XSDK_INDEX_SIZE 64")
    file(STRINGS "${SUPERLUDIST_CONFIG_PATH}" _index_size_32
         REGEX "#undef XSDK_INDEX_SIZE")
    if(_index_size_64)
      set(SUPERLUDIST_INDEX_SIZE
          64
          CACHE STRING "SuperLU DIST index size (bit width)" FORCE)
    else()
      set(SUPERLUDIST_INDEX_SIZE
          32
          CACHE STRING "SuperLU DIST index size (bit width)" FORCE)
    endif()
    mark_as_advanced(FORCE SUPERLUDIST_INDEX_SIZE)
  else()
    file(STRINGS "${SUPERLUDIST_CONFIG_PATH}" _strings_with_index_size
         REGEX "XSDK_INDEX_SIZE")
    list(GET _strings_with_index_size 0 _index_size_string)
    string(REGEX MATCHALL "[0-9][0-9]" SUPERLUDIST_INDEX_SIZE
                 "${_index_size_string}")
  endif()
  file(STRINGS "${SUPERLUDIST_CONFIG_PATH}" _strings_have_cuda
       REGEX "HAVE_CUDA")
  string(REGEX MATCH "TRUE|FALSE" _has_cuda "${_strings_have_cuda}")
  file(STRINGS "${SUPERLUDIST_CONFIG_PATH}" _strings_have_rocm REGEX "HAVE_HIP")
  string(REGEX MATCH "TRUE|FALSE" _has_rocm "${_strings_have_rocm}")
  if(_has_cuda)
    set(SUPERLUDIST_CUDA
        TRUE
        CACHE BOOL "SuperLU DIST was built with CUDA support" FORCE)
  endif()
  if(_has_rocm)
    set(SUPERLUDIST_ROCM
        TRUE
        CACHE BOOL "SuperLU DIST was built with ROCm support" FORCE)
  endif()
  unset(_has_cuda)
  unset(_has_rocm)
endif()

# find the library version file
if(NOT SUPERLUDIST_VERSION AND SUPERLUDIST_INCLUDE_DIRS)
  find_file(SUPERLUDIST_VERSION_PATH superlu_defs.h
            PATHS "${SUPERLUDIST_INCLUDE_DIRS}")

  file(STRINGS "${SUPERLUDIST_VERSION_PATH}" _version_major
       REGEX "SUPERLU_DIST_MAJOR_VERSION")
  list(GET _version_major 0 _version_string)
  string(REGEX MATCHALL "[0-9]" _version_major "${_version_string}")

  file(STRINGS "${SUPERLUDIST_VERSION_PATH}" _version_minor
       REGEX "SUPERLU_DIST_MINOR_VERSION")
  list(GET _version_minor 0 _version_string)
  string(REGEX MATCHALL "[0-9]" _version_minor "${_version_string}")

  file(STRINGS "${SUPERLUDIST_VERSION_PATH}" _version_patch
       REGEX "SUPERLU_DIST_PATCH_VERSION")
  list(GET _version_patch 0 _version_string)
  string(REGEX MATCHALL "[0-9]" _version_patch "${_version_string}")

  set(SUPERLUDIST_VERSION
      "${_version_major}.${_version_minor}.${_version_patch}")
  mark_as_advanced(FORCE SUPERLUDIST_VERSION_PATH)
endif()

# add libraries for OpenMP support
if(SUPERLUDIST_OpenMP)
  list(APPEND SUPERLUDIST_LINK_LIBRARIES OpenMP::OpenMP_CXX)
endif()
# add libraries for GPU support
if(SUPERLUDIST_CUDA)
  list(APPEND SUPERLUDIST_LINK_LIBRARIES CUDA::cudart CUDA::cublas)
endif()
if(SUPERLUDIST_ROCM)
  find_package(hipblas REQUIRED)
  find_package(rocsolver REQUIRED)
  find_package(rocblas REQUIRED)
  list(APPEND SUPERLUDIST_LINK_LIBRARIES hip::device roc::hipblas roc::rocblas
       roc::rocsolver)
endif()

# set package variables including SUPERLUDIST_FOUND
find_package_handle_standard_args(
  SUPERLUDIST
  REQUIRED_VARS SUPERLUDIST_LINK_LIBRARIES SUPERLUDIST_INCLUDE_DIRS
                SUPERLUDIST_INDEX_SIZE
  VERSION_VAR SUPERLUDIST_VERSION)

# Create target for SuperLU_DIST
if(SUPERLUDIST_FOUND)

  if(NOT TARGET SUNDIALS::SUPERLUDIST)
    add_library(SUNDIALS::SUPERLUDIST INTERFACE IMPORTED)
  endif()

  set_target_properties(
    SUNDIALS::SUPERLUDIST
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${SUPERLUDIST_INCLUDE_DIRS}"
               INTERFACE_LINK_LIBRARIES "${SUPERLUDIST_LINK_LIBRARIES}")

endif()
