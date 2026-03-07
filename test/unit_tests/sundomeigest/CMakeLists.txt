# ---------------------------------------------------------------
# Programmer(s): Mustafa Aggul @ SMU
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
# test/unit_tests/sundomeigest level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

if(SUNDIALS_ENABLE_CALIPER AND SUNDIALS_BUILD_WITH_PROFILING)
  include_directories(${caliper_INCLUDE_DIR})
  set(EXE_EXTRA_LINK_LIBS ${EXE_EXTRA_LINK_LIBS} caliper)
endif()

# Build the sundials_sundomeigest test utilities
add_library(test_sundomeigest_obj OBJECT test_sundomeigest.c
                                         test_sundomeigest.h)
if(BUILD_SHARED_LIBS)
  # need PIC when shared libs are used since the example executables will link
  # to the shared lib
  set_property(TARGET test_sundomeigest_obj PROPERTY POSITION_INDEPENDENT_CODE
                                                     TRUE)
endif()
target_link_libraries(test_sundomeigest_obj PRIVATE sundials_sundomeigestpower)

# Always add the power iteration examples
add_subdirectory(Power)

# Add the Arnoldi iteration examples if LAPACK is enabled
if(SUNDIALS_ENABLE_LAPACK)
  add_subdirectory(Arnoldi)
endif()
