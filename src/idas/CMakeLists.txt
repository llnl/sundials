# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ UMBC
#                Radu Serban, Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the IDAS library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall IDAS\n\")")

# Add variable idas_SOURCES with the sources for the IDAS library
set(idas_SOURCES
    idas.c
    idaa.c
    idas_cli.c
    idas_io.c
    idas_ic.c
    idaa_io.c
    idas_ls.c
    idas_bbdpre.c
    idas_nls.c
    idas_nls_sim.c
    idas_nls_stg.c)

# Add variable idas_HEADERS with the exported IDAS header files
set(idas_HEADERS idas.h idas.hpp idas_bbdpre.h idas_ls.h)

# Add prefix with complete path to the IDAS header files
add_prefix(${SUNDIALS_SOURCE_DIR}/include/idas/ idas_HEADERS)

# Create the library
sundials_add_library(
  sundials_idas
  SOURCES ${idas_SOURCES}
  HEADERS ${idas_HEADERS}
  INCLUDE_SUBDIR idas
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
    sundials_sunnonlinsolnewton_obj
    sundials_sunnonlinsolfixedpoint_obj
  OUTPUT_NAME sundials_idas
  VERSION ${idaslib_VERSION}
  SOVERSION ${idaslib_SOVERSION})

# Finished IDAS
message(STATUS "Added IDAS module")

# Add F2003 module if the interface is enabled
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory("fmod_int${SUNDIALS_INDEX_SIZE}")
endif()
