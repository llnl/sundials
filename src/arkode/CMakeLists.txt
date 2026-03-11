# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ UMBC, Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the ARKODE library
# ---------------------------------------------------------------

install(CODE "MESSAGE(\"\nInstall ARKODE\n\")")

# Add variable arkode_SOURCES with the sources for the ARKODE library
set(arkode_SOURCES
    arkode_adapt.c
    arkode_arkstep_io.c
    arkode_arkstep_nls.c
    arkode_arkstep.c
    arkode_bandpre.c
    arkode_bbdpre.c
    arkode_butcher_dirk.c
    arkode_butcher_erk.c
    arkode_butcher.c
    arkode_cli.c
    arkode_erkstep_io.c
    arkode_erkstep.c
    arkode_forcingstep.c
    arkode_interp.c
    arkode_io.c
    arkode_ls.c
    arkode_lsrkstep_io.c
    arkode_lsrkstep.c
    arkode_mri_tables.c
    arkode_mristep_controller.c
    arkode_mristep_io.c
    arkode_mristep_nls.c
    arkode_mristep.c
    arkode_relaxation.c
    arkode_root.c
    arkode_splittingstep_coefficients.c
    arkode_splittingstep.c
    arkode_sprkstep_io.c
    arkode_sprkstep.c
    arkode_sprk.c
    arkode_sunstepper.c
    arkode_user_controller.c
    arkode.c)

# Add variable arkode_HEADERS with the exported ARKODE header files
set(arkode_HEADERS
    arkode.h
    arkode.hpp
    arkode_arkstep.h
    arkode_arkstep.hpp
    arkode_arkstep_deprecated.h
    arkode_bandpre.h
    arkode_bbdpre.h
    arkode_butcher.h
    arkode_butcher_dirk.h
    arkode_butcher_erk.h
    arkode_erkstep.h
    arkode_erkstep.hpp
    arkode_erkstep_deprecated.h
    arkode_forcingstep.h
    arkode_forcingstep.hpp
    arkode_ls.h
    arkode_lsrkstep.h
    arkode_lsrkstep.hpp
    arkode_mristep.h
    arkode_mristep.hpp
    arkode_mristep_deprecated.h
    arkode_splittingstep.h
    arkode_splittingstep.hpp
    arkode_sprk.h
    arkode_sprkstep.h
    arkode_sprkstep.hpp
    arkode_sprkstep_deprecated.h)

# Add prefix with complete path to the ARKODE header files
add_prefix(${SUNDIALS_SOURCE_DIR}/include/arkode/ arkode_HEADERS)

# Create the sundials_arkode library
sundials_add_library(
  sundials_arkode
  SOURCES ${arkode_SOURCES}
  HEADERS ${arkode_HEADERS}
  INCLUDE_SUBDIR arkode
  LINK_LIBRARIES PUBLIC sundials_core
  OBJECT_LIBRARIES
    sundials_sunmemsys_obj
    sundials_nvecserial_obj
    sundials_nvecmanyvector_obj
    sundials_sunadaptcontrollersoderlind_obj
    sundials_sunadaptcontrollerimexgus_obj
    sundials_sunadaptcontrollermrihtol_obj
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
    sundials_adjointcheckpointscheme_fixed_obj
  OUTPUT_NAME sundials_arkode
  VERSION ${arkodelib_VERSION}
  SOVERSION ${arkodelib_SOVERSION})

# Finished ARKODE
message(STATUS "Added ARKODE module")

# Add F2003 module if the interface is enabled
if(SUNDIALS_ENABLE_FORTRAN)
  add_subdirectory("fmod_int${SUNDIALS_INDEX_SIZE}")
endif()

# Add ARKODE XBraid interface
if(SUNDIALS_ENABLE_XBRAID)
  add_subdirectory(xbraid)
endif()
