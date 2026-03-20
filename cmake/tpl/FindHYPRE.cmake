# ---------------------------------------------------------------
# Programmer(s): Eddy Banks, Slaven Peles, Cody J. Balos, and
#                Jean Sexton @ LLNL
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
# HYPRE find module that creates an imported target for HYPRE.
# The target is SUNDIALS::HYPRE.
#
# The variable HYPRE_LIBRARY_DIR can be used to control
# where the module looks for the library.
#
# The variable HYPRE_INCLUDE_DIR can be used to set the
# include path for the library.
#
# This module also defines variables, but it is best to use
# the defined target to ensure includes and compile/link
# options are correctly passed to consumers.
#
#   HYPRE_FOUND       - system has HYPRE library
#   HYPRE_DIR
#   HYPRE_LIBRARY     - the HYPRE library
#   HYPRE_INCLUDE_DIR - the HYPRE include path
#   HYPRE_LIBRARIES   - all of the libraries needed for HYPRE
# ---------------------------------------------------------------

# Find include dir
find_path(
  temp_HYPRE_INCLUDE_DIR
  NAMES HYPRE.h hypre.h
  HINTS "${HYPRE_DIR}" "${HYPRE_DIR}/include" "${HYPRE_INCLUDE_DIR}")
if(temp_HYPRE_INCLUDE_DIR)
  set(HYPRE_INCLUDE_DIR
      "${temp_HYPRE_INCLUDE_DIR}"
      CACHE PATH "" FORCE)
endif()
unset(temp_HYPRE_INCLUDE_DIR CACHE)

if(HYPRE_LIBRARY)
  # We have (or were given) HYPRE_LIBRARY - get path to use for any related libs
  get_filename_component(HYPRE_LIBRARY_DIR ${HYPRE_LIBRARY} PATH)

  # force CACHE update to show user DIR that will be used
  set(HYPRE_LIBRARY_DIR
      ${HYPRE_LIBRARY_DIR}
      CACHE PATH "" FORCE)
else()
  # find library with user provided directory path
  set(HYPRE_LIBRARY_NAMES hypre HYPRE)
  find_library(
    HYPRE_LIBRARY
    NAMES ${HYPRE_LIBRARY_NAMES}
    HINTS "${HYPRE_DIR}" "${HYPRE_DIR}/lib" "${HYPRE_DIR}/lib64"
          "${HYPRE_LIBRARY_DIR}")
endif()
mark_as_advanced(HYPRE_LIBRARY)

list(FIND HYPRE_LIBRARIES ${HYPRE_LIBRARY} _idx)
if(_idx EQUAL -1)
  set(HYPRE_LIBRARIES
      "${HYPRE_LIBRARY};${HYPRE_LIBRARIES}"
      CACHE STRING "" FORCE)
endif()

# Extract version parts from the version information
if(HYPRE_INCLUDE_DIR)

  # HYPRE_config.h file added in at least v2.11.0 (possibly sooner)
  find_file(HYPRE_CONFIG_FILE HYPRE_config.h PATHS "${HYPRE_INCLUDE_DIR}")

  file(STRINGS "${HYPRE_CONFIG_FILE}" _hypre_release_version
       REGEX "HYPRE_RELEASE_VERSION")
  string(REGEX MATCH "[0-9]+\.[0-9]+\.[0-9]+" HYPRE_VERSION
               "${_hypre_release_version}")
  string(REGEX MATCHALL "[0-9]+" _hypre_version_numbers "${HYPRE_VERSION}")
  list(GET _hypre_version_numbers 0 HYPRE_VERSION_MAJOR)
  list(GET _hypre_version_numbers 1 HYPRE_VERSION_MINOR)
  list(GET _hypre_version_numbers 2 HYPRE_VERSION_PATCH)

endif()

# set a more informative error message in case the library was not found
set(HYPRE_NOT_FOUND_MESSAGE
    "\
************************************************************************\n\
ERROR: Could not find HYPRE. Please check the variables:\n\
       HYPRE_INCLUDE_DIR and HYPRE_LIBRARY_DIR\n\
************************************************************************")

# set package variables including HYPRE_FOUND
find_package_handle_standard_args(
  HYPRE
  REQUIRED_VARS HYPRE_LIBRARY HYPRE_LIBRARIES HYPRE_INCLUDE_DIR
  VERSION_VAR HYPRE_VERSION
  FAIL_MESSAGE "${HYPRE_NOT_FOUND_MESSAGE}")

# Create target for HYPRE
if(HYPRE_FOUND)

  if(NOT TARGET SUNDIALS::HYPRE)
    add_library(SUNDIALS::HYPRE UNKNOWN IMPORTED)
  endif()

  set_target_properties(
    SUNDIALS::HYPRE
    PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${HYPRE_INCLUDE_DIR}"
               INTERFACE_LINK_LIBRARIES "${HYPRE_LIBRARIES}"
               IMPORTED_LOCATION "${HYPRE_LIBRARY}")

endif()
