# -----------------------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------------------
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
# -----------------------------------------------------------------------------
# Module to find and setup OpenMP.
#
# Creates the variables:
#   OPENMP_FOUND - was OpenMP found
#   OPENMP45_FOUND - was OpenMP v4.5 or greater found
#   OPENMP_SUPPORTS_DEVICE_OFFLOADING - is device offloading supported
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 1: Include guard
# -----------------------------------------------------------------------------

include_guard(GLOBAL)

# -----------------------------------------------------------------------------
# Section 2: Check to make sure options are compatible
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Section 3: Find the TPL
# -----------------------------------------------------------------------------

set(OPENMP45_FOUND FALSE)
set(OPENMP_SUPPORTS_DEVICE_OFFLOADING FALSE)
find_package(OpenMP REQUIRED)

# -----------------------------------------------------------------------------
# Section 4: Test the TPL
# -----------------------------------------------------------------------------

# Work around a bug in setting OpenMP version variables in CMake >= 3.9. The
# OpenMP version information is not stored in cache variables and is not set on
# repeated calls to find OpenMP (i.e., when using ccmake). To ensure these
# variables exist store copies of the values.
set(OpenMP_C_VERSION
    "${OpenMP_C_VERSION}"
    CACHE INTERNAL "" FORCE)
set(OpenMP_CXX_VERSION
    "${OpenMP_CXX_VERSION}"
    CACHE INTERNAL "" FORCE)
set(OpenMP_Fortran_VERSION
    "${OpenMP_Fortran_VERSION}"
    CACHE INTERNAL "" FORCE)

# Check for OpenMP offloading support
if(OPENMP_FOUND AND (SUNDIALS_ENABLE_OPENMP_DEVICE OR SUPERLUDIST_OpenMP))

  if(SUNDIALS_ENABLE_OPENMP_DEVICE_CHECKS)

    # Check the OpenMP version
    message(STATUS "Checking whether OpenMP supports device offloading")

    if((OpenMP_C_VERSION VERSION_EQUAL 4.5) OR (OpenMP_C_VERSION VERSION_GREATER
                                                4.5))
      message(
        STATUS "Checking whether OpenMP supports device offloading -- yes")
      set(OPENMP45_FOUND TRUE)
      set(OPENMP_SUPPORTS_DEVICE_OFFLOADING TRUE)
    else()
      message(STATUS "Checking whether OpenMP supports device offloading -- no")
      set(OPENMP45_FOUND FALSE)
      set(OPENMP_SUPPORTS_DEVICE_OFFLOADING FALSE)
      message(
        FATAL_ERROR
          "The found OpenMP version does not support device offloading.")
    endif()

  else()

    # The user has asked for checks to be skipped, assume offloading is
    # supported
    set(OPENMP45_FOUND TRUE)
    set(OPENMP_SUPPORTS_DEVICE_OFFLOADING TRUE)
    message(STATUS "Skipped OpenMP checks.")
    message(
      WARNING
        "SUNDIALS OpenMP functionality dependent on OpenMP 4.5+ is not guaranteed."
    )

  endif()

endif()
