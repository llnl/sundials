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

# Create a library out of the generic sundials modules
sundials_add_library(
  sundials_adjointcheckpointscheme_fixed
  SOURCES sunadjointcheckpointscheme_fixed.c
  HEADERS
    ${SUNDIALS_SOURCE_DIR}/include/sundials/sundials_adjointstepper.h
    ${SUNDIALS_SOURCE_DIR}/include/sundials/sundials_adjointcheckpointscheme.h
    ${SUNDIALS_SOURCE_DIR}/include/sunadjointcheckpointscheme/sunadjointcheckpointscheme_fixed.h
  LINK_LIBRARIES PUBLIC sundials_core
  INCLUDE_SUBDIR sunadjointcheckpointscheme
  OBJECT_LIB_ONLY)

# Add F2003 module if the interface is enabled
if(SUNDIALS_ENABLE_FORTRAN)
  add_subdirectory("fmod_int${SUNDIALS_INDEX_SIZE}")
endif()
