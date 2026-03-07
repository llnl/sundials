# ---------------------------------------------------------------
# Programmer(s): David J. Gardner, Slaven Peles, and
#                Cody Balos @ LLNL
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
# test/unit_tests/sunmatrix level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

if(SUNDIALS_ENABLE_CALIPER AND SUNDIALS_BUILD_WITH_PROFILING)
  include_directories(${caliper_INCLUDE_DIR})
  set(EXE_EXTRA_LINK_LIBS ${EXE_EXTRA_LINK_LIBS} caliper)
endif()

# Always add the serial sunmatrix dense/band/sparse examples
add_subdirectory(dense)
add_subdirectory(band)
add_subdirectory(sparse)

# Build the sunmatrix test utilities
add_library(test_sunmatrix_obj OBJECT test_sunmatrix.c test_sunmatrix.h)
if(BUILD_SHARED_LIBS)
  # need PIC when shared libs are used since the example executables will link
  # to the shared lib
  set_property(TARGET test_sunmatrix_obj PROPERTY POSITION_INDEPENDENT_CODE
                                                  TRUE)
endif()
target_link_libraries(test_sunmatrix_obj PRIVATE sundials_sunmatrixdense)

if(BUILD_SUNMATRIX_CUSPARSE)
  if(SUNDIALS_INDEX_SIZE MATCHES "32")
    add_subdirectory(cusparse)
  endif()
endif()

if(BUILD_SUNMATRIX_GINKGO)
  add_subdirectory(ginkgo)
endif()

if(BUILD_SUNMATRIX_KOKKOSDENSE)
  add_subdirectory(kokkos)
endif()

if(BUILD_SUNMATRIX_MAGMADENSE)
  add_subdirectory(magmadense)
endif()

if(BUILD_SUNMATRIX_ONEMKLDENSE)
  add_subdirectory(onemkldense)
endif()

if(BUILD_SUNMATRIX_SLUNRLOC)
  add_subdirectory(slunrloc)
endif()
