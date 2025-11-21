# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ UMBC
#                Radu Serban, Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025, Lawrence Livermore National Security,
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
# CMakeLists.txt file for the KINSOL library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall KINSOL\n\")")

# Add variable kinsol_SOURCES with the sources for the KINSOL library
set(kinsol_SOURCES
    kinsol_aa.c
    kinsol_bbdpre.c
    kinsol_cli.c
    kinsol_io.c
    kinsol_ls.c
    kinsol_orth.c
    kinsol.c)

# Add variable kinsol_HEADERS with the exported KINSOL header files
set(kinsol_HEADERS kinsol.h kinsol.hpp kinsol_bbdpre.h kinsol_ls.h)

# Add prefix with complete path to the KINSOL header files
add_prefix(${SUNDIALS_SOURCE_DIR}/include/kinsol/ kinsol_HEADERS)

# Create the library
sundials_add_library(
  sundials_kinsol
  SOURCES ${kinsol_SOURCES}
  HEADERS ${kinsol_HEADERS}
  INCLUDE_SUBDIR kinsol
  LINK_LIBRARIES PUBLIC sundials_core
  OBJECT_LIBRARIES
    sundials_sunmemsys_obj
    sundials_nvecserial_obj
    sundials_sunmatrixband_obj
    sundials_sunmatrixdense_obj
    sundials_sunmatrixsparse_obj
    sundials_sunlinsolband_obj
    sundials_sunlinsoldense_obj
    sundials_sunlinsolspbcgs_obj
    sundials_sunlinsolspfgmr_obj
    sundials_sunlinsolspgmr_obj
    sundials_sunlinsolsptfqmr_obj
    sundials_sunlinsolpcg_obj
  OUTPUT_NAME sundials_kinsol
  VERSION ${kinsollib_VERSION}
  SOVERSION ${kinsollib_SOVERSION})

# Finished KINSOL
message(STATUS "Added KINSOL module")

# Add F2003 module if the interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory("fmod_int${SUNDIALS_INDEX_SIZE}")
endif()
