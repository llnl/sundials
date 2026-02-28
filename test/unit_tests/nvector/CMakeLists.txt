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
# test/unit_tests/nvector level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

if(ENABLE_CALIPER AND SUNDIALS_BUILD_WITH_PROFILING)
  include_directories(${caliper_INCLUDE_DIR})
  set(EXE_EXTRA_LINK_LIBS ${EXE_EXTRA_LINK_LIBS} caliper)
endif()

# Always add the nvector serial examples
add_subdirectory(serial)

# Build the nvector test utilities
add_library(test_nvector_obj OBJECT test_nvector.c)
if(BUILD_SHARED_LIBS)
  # need PIC when shared libs are used since the example executables will link
  # to the shared lib
  set_property(TARGET test_nvector_obj PROPERTY POSITION_INDEPENDENT_CODE TRUE)
endif()
target_link_libraries(test_nvector_obj PRIVATE sundials_nvecserial)

if(ENABLE_MPI AND MPI_C_FOUND)
  add_subdirectory(parallel)
  if(BUILD_NVECTOR_MPIMANYVECTOR)
    add_subdirectory(mpimanyvector)
  endif()
  add_subdirectory(mpiplusx)
  # Build the mpi nvector test utilities
  add_library(test_nvectormpi_obj OBJECT test_mpinvector.c)
  if(BUILD_SHARED_LIBS)
    set_property(TARGET test_nvectormpi_obj PROPERTY POSITION_INDEPENDENT_CODE
                                                     TRUE)
  endif()
  target_link_libraries(test_nvectormpi_obj PRIVATE MPI::MPI_C
                                                    sundials_nvecparallel)
endif()
if(BUILD_NVECTOR_MANYVECTOR)
  add_subdirectory(manyvector)
endif()

if(BUILD_NVECTOR_PARHYP)
  add_subdirectory(parhyp)
endif()

if(BUILD_NVECTOR_PTHREADS)
  add_subdirectory(pthreads)
endif()

if(BUILD_NVECTOR_OPENMP)
  add_subdirectory(C_openmp)
endif()

if(BUILD_NVECTOR_OPENMPDEV)
  add_subdirectory(openmpdev)
endif()

if(BUILD_NVECTOR_PETSC)
  add_subdirectory(petsc)
endif()

if(BUILD_NVECTOR_CUDA)
  add_subdirectory(cuda)
  if(ENABLE_MPI AND MPI_C_FOUND)
    add_subdirectory(mpicuda)
  endif()
endif()

if(BUILD_NVECTOR_HIP)
  add_subdirectory(hip)
endif()

if(BUILD_NVECTOR_SYCL)
  add_subdirectory(sycl)
endif()

if(BUILD_NVECTOR_RAJA)
  add_subdirectory(raja)
  if(ENABLE_MPI AND MPI_C_FOUND)
    add_subdirectory(mpiraja)
  endif()
endif()

if(BUILD_NVECTOR_TRILINOS)
  add_subdirectory(trilinos)
endif()

if(BUILD_NVECTOR_KOKKOS)
  add_subdirectory(kokkos)
endif()
