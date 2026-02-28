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

find_package(PkgConfig REQUIRED)

set(_petsc_prefixes "")
if(DEFINED PETSC_DIR)
  list(APPEND _petsc_prefixes "${PETSC_DIR}")
  if(DEFINED PETSC_ARCH)
    list(APPEND _petsc_prefixes "${PETSC_DIR}/${PETSC_ARCH}")
  endif()
endif()
list(APPEND CMAKE_PREFIX_PATH ${_petsc_prefixes})
unset(_petsc_prefixes)

# Use pkg-config to find PETSC
set(PKG_CONFIG_USE_CMAKE_PREFIX_PATH "YES")
set(_pkg_version_spec "")
if(DEFINED PETSC_FIND_VERSION)
  if(PETSC_FIND_VERSION_EXACT)
    set(_pkg_version_spec "=${PETSC_FIND_VERSION}")
  else()
    set(_pkg_version_spec ">=${PETSC_FIND_VERSION}")
  endif()
endif()
pkg_check_modules(PKG_PETSC "PETSc${_pkg_version_spec}")
unset(_pkg_version_spec)

# Find the PETSC libraries
set(_petsc_libs)
foreach(_next_lib IN LISTS PKG_PETSC_LIBRARIES)
  find_library(
    _petsc_lib_${_next_lib}
    NAMES ${_next_lib}
    HINTS ${PKG_PETSC_LIBRARY_DIRS})
  if(_petsc_lib_${_next_lib})
    list(APPEND _petsc_libs "${_petsc_lib_${_next_lib}}")
  endif()
endforeach()

# libm is always required
list(APPEND _petsc_libs "${SUNDIALS_MATH_LIBRARY}")

# Substitute MPI target if PETSC is built with MPI
foreach(_next_lib IN LISTS PKG_PETSC_STATIC_LIBRARIES)
  if(_next_lib MATCHES "mpi")
    list(APPEND _petsc_libs "MPI::MPI_C")
  endif()
  if(_next_lib MATCHES "kokkoskernels")
    if(NOT TARGET Kokkos::kokkoskernels)
      find_package(KokkosKernels REQUIRED HINTS "${KokkosKernels_DIR}"
                   "${PKG_PETSC_LIBRARY_DIRS}" NO_DEFAULT_PATH)
    endif()
    list(APPEND _petsc_libs "Kokkos::kokkoskernels")
  endif()
  if(_next_lib MATCHES "kokkos")
    if(NOT TARGET Kokkos::kokkos)
      find_package(Kokkos REQUIRED HINTS "${Kokkos_DIR}"
                   "${PKG_PETSC_LIBRARY_DIRS}" NO_DEFAULT_PATH)
    endif()
    list(APPEND _petsc_libs "Kokkos::kokkos")
  endif()
endforeach()
list(REMOVE_DUPLICATES _petsc_libs)

# Set result variables
set(PETSC_LIBRARIES "${_petsc_libs}")
unset(_petsc_libs)
set(PETSC_FOUND ${PKG_PETSC_FOUND})
set(PETSC_INCLUDE_DIRS ${PKG_PETSC_INCLUDE_DIRS})

# Extract version parts from the version information
if(PKG_PETSC_VERSION)
  set(_petsc_versions "")
  string(REGEX MATCHALL "[0-9]+" _petsc_versions ${PKG_PETSC_VERSION})
  list(GET _petsc_versions 0 _petsc_version_major)
  list(GET _petsc_versions 1 _petsc_version_minor)
  list(GET _petsc_versions 2 _petsc_version_patch)

  set(PETSC_VERSION
      ${PKG_PETSC_VERSION}
      CACHE STRING "Full version of PETSC")
  set(PETSC_VERSION_MAJOR
      ${_petsc_version_major}
      CACHE INTERNAL "Major version of PETSC")
  set(PETSC_VERSION_MINOR
      ${_petsc_version_minor}
      CACHE INTERNAL "Minor version of PETSC")
  set(PETSC_VERSION_PATCH
      ${_petsc_version_patch}
      CACHE INTERNAL "Patch version of PETSC")

  unset(_petsc_versions)
  unset(_petsc_version_major)
  unset(_petsc_version_minor)
  unset(_petsc_version_patch)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  PETSC
  REQUIRED_VARS PETSC_FOUND PETSC_INCLUDE_DIRS PETSC_LIBRARIES
  VERSION_VAR PETSC_VERSION)

if(NOT TARGET SUNDIALS::PETSC)
  add_library(SUNDIALS::PETSC INTERFACE IMPORTED)
  set_target_properties(
    SUNDIALS::PETSC
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${PETSC_INCLUDE_DIRS}"
               INTERFACE_LINK_LIBRARIES "${PETSC_LIBRARIES}")
endif()

mark_as_advanced(PETSC_INCLUDE_DIRS PETSC_LIBRARIES PETSC_VERSION_MAJOR
                 PETSC_VERSION_MINOR PETSC_VERSION_PATCH PETSC_VERSION)
