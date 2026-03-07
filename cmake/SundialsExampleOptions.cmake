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
if(SUNDIALS_ENABLE_TRILINOS
   OR SUNDIALS_ENABLE_SUPERLUDIST
   OR SUNDIALS_ENABLE_XBRAID
   OR SUNDIALS_ENABLE_HIP
   OR SUNDIALS_ENABLE_MAGMA
   OR SUNDIALS_ENABLE_SYCL
   OR SUNDIALS_ENABLE_ONEMKL
   OR SUNDIALS_ENABLE_RAJA
   OR SUNDIALS_ENABLE_GINKGO
   OR SUNDIALS_ENABLE_KOKKOS)
  sundials_option(EXAMPLES_ENABLE_CXX BOOL "Build SUNDIALS C++ examples" ON)
else()
  sundials_option(EXAMPLES_ENABLE_CXX BOOL "Build SUNDIALS C++ examples" OFF)
endif()

# -----------------------------------------------------------------------------
# Options for Fortran Examples
# -----------------------------------------------------------------------------

if(SUNDIALS_ENABLE_FORTRAN)
  sundials_option(EXAMPLES_ENABLE_F2003 BOOL "Build SUNDIALS Fortran examples"
                  ON DEPENDS_ON SUNDIALS_ENABLE_FORTRAN)
else()
  sundials_option(EXAMPLES_ENABLE_F2003 BOOL "Build SUNDIALS Fortran examples"
                  OFF DEPENDS_ON SUNDIALS_ENABLE_FORTRAN)
endif()

# -----------------------------------------------------------------------------
# Options for CUDA Examples
# -----------------------------------------------------------------------------

sundials_option(EXAMPLES_ENABLE_CUDA BOOL "Build SUNDIALS CUDA examples" ON
                DEPENDS_ON SUNDIALS_ENABLE_CUDA)

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
