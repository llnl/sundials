# ---------------------------------------------------------------
# Programmer(s): David J. Gardner and Slaven Peles @ LLNL
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
# test/unit_tests/sunlinsol level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

if(ENABLE_CALIPER AND SUNDIALS_BUILD_WITH_PROFILING)
  include_directories(${caliper_INCLUDE_DIR})
  set(EXE_EXTRA_LINK_LIBS ${EXE_EXTRA_LINK_LIBS} caliper)
endif()

# Always add the serial sunlinearsolver dense and band examples
add_subdirectory(band)
add_subdirectory(dense)

# Always add serial sunlinearsolver iterative examples
add_subdirectory(spgmr/serial)
add_subdirectory(spfgmr/serial)
add_subdirectory(spbcgs/serial)
add_subdirectory(sptfqmr/serial)
add_subdirectory(pcg/serial)

# Build the sunlinsol test utilities
add_library(test_sunlinsol_obj OBJECT test_sunlinsol.c test_sunlinsol.h)
if(BUILD_SHARED_LIBS)
  # need PIC when shared libs are used since the example executables will link
  # to the shared lib
  set_property(TARGET test_sunlinsol_obj PROPERTY POSITION_INDEPENDENT_CODE
                                                  TRUE)
endif()
target_link_libraries(test_sunlinsol_obj PRIVATE sundials_sunlinsoldense)

if(ENABLE_MPI AND MPI_C_FOUND)
  add_subdirectory(pcg/parallel)
  add_subdirectory(spgmr/parallel)
  add_subdirectory(spfgmr/parallel)
  add_subdirectory(spbcgs/parallel)
  add_subdirectory(sptfqmr/parallel)
endif()

if(BUILD_SUNMATRIX_CUSPARSE)
  if(SUNDIALS_INDEX_SIZE MATCHES "32")
    add_subdirectory(cusolversp)
  endif()
endif()

if(BUILD_SUNLINSOL_GINKGO)
  add_subdirectory(ginkgo)
endif()

if(BUILD_SUNLINSOL_KLU)
  add_subdirectory(klu)
endif()

if(BUILD_SUNLINSOL_KOKKOSDENSE)
  add_subdirectory(kokkos)
endif()

if(BUILD_SUNLINSOL_LAPACKBAND)
  add_subdirectory(lapackband)
endif()

if(BUILD_SUNLINSOL_LAPACKDENSE)
  add_subdirectory(lapackdense)
endif()

if(BUILD_SUNLINSOL_MAGMADENSE)
  add_subdirectory(magmadense)
endif()

if(BUILD_SUNLINSOL_ONEMKLDENSE)
  add_subdirectory(onemkldense)
endif()

if(BUILD_SUNLINSOL_SUPERLUMT)
  add_subdirectory(superlumt)
endif()

if(BUILD_SUNLINSOL_SUPERLUDIST)
  add_subdirectory(superludist)
endif()
