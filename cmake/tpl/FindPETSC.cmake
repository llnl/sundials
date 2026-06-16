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
    PATHS ${PKG_PETSC_LIBRARY_DIRS})
  if(_petsc_lib_${_next_lib})
    list(APPEND _petsc_libs "${_petsc_lib_${_next_lib}}")
  endif()
endforeach()

# libm is always required
list(APPEND _petsc_libs "${SUNDIALS_MATH_LIBRARY}")

# Is libpetsc static?
set(_petsc_static FALSE)
if(_petsc_lib_petsc AND _petsc_lib_petsc MATCHES "\\.(a|lib)$")
  set(_petsc_static TRUE)
endif()

set(_petsc_link_options "")
if(_petsc_static)
  # Names already resolved from the public Libs (skip to avoid double-link)
  set(_petsc_public ${PKG_PETSC_LIBRARIES})
  foreach(_flag IN LISTS PKG_PETSC_STATIC_LDFLAGS)
    if(_flag MATCHES "^-l(.+)$")
      if(CMAKE_MATCH_1 IN_LIST _petsc_public)
        continue()
      endif()
    endif()
    if(_flag MATCHES "^-L")
      list(APPEND _petsc_link_options "${_flag}")
    elseif(_flag MATCHES "^-l" OR _flag MATCHES "^/")
      list(APPEND _petsc_libs "${_flag}")
    else()
      list(APPEND _petsc_link_options "${_flag}")
    endif()
  endforeach()
  list(APPEND _petsc_link_options ${PKG_PETSC_STATIC_LDFLAGS_OTHER})
  unset(_petsc_public)
endif()

# Set result variables
set(PETSC_LIBRARIES "${_petsc_libs}")
unset(_petsc_libs)
set(PETSC_LINK_OPTIONS "${_petsc_link_options}")
unset(_petsc_link_options)
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
  if(PETSC_LINK_OPTIONS)
    set_target_properties(SUNDIALS::PETSC PROPERTIES
                          INTERFACE_LINK_OPTIONS "${PETSC_LINK_OPTIONS}")
  endif()
endif()

mark_as_advanced(PETSC_INCLUDE_DIRS PETSC_LIBRARIES PETSC_VERSION_MAJOR
                 PETSC_VERSION_MINOR PETSC_VERSION_PATCH PETSC_VERSION)
