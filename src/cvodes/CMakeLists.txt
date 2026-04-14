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
# CMakeLists.txt file for the CVODES library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall CVODES\n\")")

# Add variable cvodes_SOURCES with the sources for the CVODES library
set(cvodes_SOURCES
    cvodea.c
    cvodea_io.c
    cvodes.c
    cvodes_bandpre.c
    cvodes_bbdpre.c
    cvodes_cli.c
    cvodes_diag.c
    cvodes_io.c
    cvodes_ls.c
    cvodes_nls.c
    cvodes_nls_sim.c
    cvodes_nls_stg.c
    cvodes_nls_stg1.c
    cvodes_proj.c
    cvodes_resize.c)

# Add variable cvodes_HEADERS with the exported CVODES header files
set(cvodes_HEADERS
    cvodes.h
    cvodes.hpp
    cvodes_bandpre.h
    cvodes_bbdpre.h
    cvodes_diag.h
    cvodes_ls.h
    cvodes_proj.h)

# Add prefix with complete path to the CVODES header files
add_prefix(${SUNDIALS_SOURCE_DIR}/include/cvodes/ cvodes_HEADERS)

# Create the library
sundials_add_library(
  sundials_cvodes
  SOURCES ${cvodes_SOURCES}
  HEADERS ${cvodes_HEADERS}
  INCLUDE_SUBDIR cvodes
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
  OUTPUT_NAME sundials_cvodes
  VERSION ${cvodeslib_VERSION}
  SOVERSION ${cvodeslib_SOVERSION})

# Finished CVODES
message(STATUS "Added CVODES module")

# Add F2003 module if the interface is enabled
if(SUNDIALS_ENABLE_FORTRAN)
  add_subdirectory("fmod_int${SUNDIALS_INDEX_SIZE}")
endif()
