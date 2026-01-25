# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# Options for SUNDIALS examples.
# ---------------------------------------------------------------

# -----------------------------------------------------------------------------
# Options for C/C++ examples
# -----------------------------------------------------------------------------

sundials_option(EXAMPLES_ENABLE_C BOOL "Build SUNDIALS C examples" ON)

# Some TPLs only have C++ examples. Default the C++ examples to ON if any of
# these are enabled on the initial configuration pass.
if(ENABLE_TRILINOS
   OR ENABLE_SUPERLUDIST
   OR ENABLE_XBRAID
   OR ENABLE_HIP
   OR ENABLE_MAGMA
   OR ENABLE_SYCL
   OR ENABLE_ONEMKL
   OR ENABLE_RAJA
   OR ENABLE_GINKGO
   OR ENABLE_KOKKOS)
  sundials_option(EXAMPLES_ENABLE_CXX BOOL "Build SUNDIALS C++ examples" ON)
else()
  sundials_option(EXAMPLES_ENABLE_CXX BOOL "Build SUNDIALS C++ examples" OFF)
endif()

# -----------------------------------------------------------------------------
# Options for Fortran Examples
# -----------------------------------------------------------------------------

# F2003 examples (on by default) are an option only if the Fortran 2003
# interface is enabled.
set(DOCSTR "Build SUNDIALS Fortran 2003 examples")
if(BUILD_FORTRAN_MODULE_INTERFACE)
  set(EXAMPLES_ENABLE_F2003
      ON
      CACHE BOOL "${DOCSTR}")

  # Fortran 2003 examples only support double precision
  if(EXAMPLES_ENABLE_F2003 AND (NOT (SUNDIALS_PRECISION MATCHES "DOUBLE")))
    message(
      WARNING
        "F2003 examples are not compatible with ${SUNDIALS_PRECISION} precision. "
        "Setting EXAMPLES_ENABLE_F2003 to OFF.")
    set(EXAMPLES_ENABLE_F2003
        OFF
        CACHE BOOL "${DOCSTR}" FORCE)
  endif()
else()

  # set back to OFF (in case it was ON)
  if(EXAMPLES_ENABLE_F2003)
    message(
      WARNING
        "EXAMPLES_ENABLE_F2003 is ON but BUILD_FORTRAN_MODULE_INTERFACE is OFF. "
        "Setting EXAMPLES_ENABLE_F2003 to OFF.")
    set(EXAMPLES_ENABLE_F2003
        OFF
        CACHE BOOL "${DOCSTR}" FORCE)
  endif()

endif()

# -----------------------------------------------------------------------------
# Options for CUDA Examples
# -----------------------------------------------------------------------------

sundials_option(EXAMPLES_ENABLE_CUDA BOOL "Build SUNDIALS CUDA examples" ON
                DEPENDS_ON ENABLE_CUDA)

# -----------------------------------------------------------------------------
# Options for installing examples
# -----------------------------------------------------------------------------

# Enable installing examples by default
sundials_option(EXAMPLES_INSTALL BOOL "Install SUNDIALS examples" ON)

sundials_option(
  EXAMPLES_INSTALL_PATH PATH "Output directory for installing example files"
  "${CMAKE_INSTALL_PREFIX}/examples")

# If examples are to be exported, check where we should install them.
if(EXAMPLES_INSTALL AND NOT EXAMPLES_INSTALL_PATH)
  message(
    WARNING "The example installation path is empty. Example installation "
            "path was reset to its default value")
  set(EXAMPLES_INSTALL_PATH
      "${CMAKE_INSTALL_PREFIX}/examples"
      CACHE STRING "Output directory for installing example files" FORCE)
endif()

# -----------------------------------------------------------------------------
# Internal variables.
# -----------------------------------------------------------------------------

if(EXAMPLES_ENABLE_C
   OR EXAMPLES_ENABLE_CXX
   OR EXAMPLES_ENABLE_CUDA
   OR EXAMPLES_ENABLE_F2003)
  set(_BUILD_EXAMPLES
      TRUE
      CACHE INTERNAL "")
else()
  set(_BUILD_EXAMPLES
      FALSE
      CACHE INTERNAL "")
endif()
