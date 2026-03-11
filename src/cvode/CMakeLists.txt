# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ UMBC
#                Radu Serban, and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the CVODE library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall CVODE\n\")")

# Add variable cvode_SOURCES with the sources for the CVODE library
set(cvode_SOURCES
    cvode.c
    cvode_bandpre.c
    cvode_bbdpre.c
    cvode_cli.c
    cvode_diag.c
    cvode_io.c
    cvode_ls.c
    cvode_nls.c
    cvode_proj.c
    cvode_resize.c)

# Add variable cvode_HEADERS with the exported CVODE header files
set(cvode_HEADERS
    cvode.h
    cvode.hpp
    cvode_bandpre.h
    cvode_bbdpre.h
    cvode_diag.h
    cvode_ls.h
    cvode_proj.h)

# Add prefix with complete path to the CVODE header files
add_prefix(${SUNDIALS_SOURCE_DIR}/include/cvode/ cvode_HEADERS)

# Build fused kernel libraries
if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)

  if(BUILD_NVECTOR_CUDA)

    set_source_files_properties(cvode_fused_gpu.cpp PROPERTIES LANGUAGE CUDA)

    sundials_add_library(
      sundials_cvode_fused_cuda
      SOURCES cvode_fused_gpu.cpp
      COMPILE_DEFINITIONS PRIVATE USE_CUDA
      LINK_LIBRARIES PUBLIC sundials_core PRIVATE sundials_nveccuda
      OUTPUT_NAME sundials_cvode_fused_cuda
      VERSION ${cvodelib_VERSION}
      SOVERSION ${cvodelib_SOVERSION})
  endif()

  if(BUILD_NVECTOR_HIP)
    sundials_add_library(
      sundials_cvode_fused_hip
      SOURCES cvode_fused_gpu.cpp
      COMPILE_DEFINITIONS PRIVATE USE_HIP
      LINK_LIBRARIES PUBLIC sundials_core PRIVATE sundials_nvechip
      OUTPUT_NAME sundials_cvode_fused_hip
      VERSION ${cvodelib_VERSION}
      SOVERSION ${cvodelib_SOVERSION})
  endif()

  sundials_add_library(
    sundials_cvode_fused_stubs
    SOURCES cvode_fused_stubs.c
    LINK_LIBRARIES PUBLIC sundials_core
    OUTPUT_NAME sundials_cvode_fused_stubs
    VERSION ${cvodelib_VERSION}
    SOVERSION ${cvodelib_SOVERSION})

  set(_fused_link_lib sundials_cvode_fused_stubs)
endif()

# Create the library
sundials_add_library(
  sundials_cvode
  SOURCES ${cvode_SOURCES}
  HEADERS ${cvode_HEADERS}
  INCLUDE_SUBDIR cvode
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
  LINK_LIBRARIES # Link to stubs so examples work.
                 PRIVATE ${_fused_link_lib}
  OUTPUT_NAME sundials_cvode
  VERSION ${cvodelib_VERSION}
  SOVERSION ${cvodelib_SOVERSION})

# Finished CVODE
message(STATUS "Added CVODE module")

# Add F2003 module if the interface is enabled
if(SUNDIALS_ENABLE_FORTRAN)
  add_subdirectory("fmod_int${SUNDIALS_INDEX_SIZE}")
endif()
