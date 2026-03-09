# ---------------------------------------------------------------
# Programmer(s): Radu Serban and Cody J. Balos @ LLNL
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
# CMakeLists.txt file for the generic SUNDIALS modules
# ---------------------------------------------------------------

# From here we only install the generic SUNDIALS headers. The implementations
# themselves are incorporated in the individual SUNDIALS solver libraries.

install(CODE "MESSAGE(\"\nInstall shared components\n\")")

# Add variable sundials_HEADERS with the exported SUNDIALS header files
set(sundials_HEADERS
    sundials_adaptcontroller.h
    sundials_adaptcontroller.hpp
    sundials_adjointcheckpointscheme.h
    sundials_adjointcheckpointscheme.hpp
    sundials_adjointstepper.h
    sundials_adjointstepper.hpp
    sundials_band.h
    sundials_base.hpp
    sundials_classview.hpp
    sundials_context.h
    sundials_context.hpp
    sundials_convertibleto.hpp
    sundials_core.h
    sundials_core.hpp
    sundials_dense.h
    sundials_direct.h
    sundials_domeigestimator.h
    sundials_domeigestimator.hpp
    sundials_errors.h
    sundials_futils.h
    sundials_iterative.h
    sundials_linearsolver.h
    sundials_linearsolver.hpp
    sundials_logger.h
    sundials_logger.hpp
    sundials_math.h
    sundials_matrix.h
    sundials_matrix.hpp
    sundials_memory.h
    sundials_memory.hpp
    sundials_mpi_types.h
    sundials_nonlinearsolver.h
    sundials_nonlinearsolver.hpp
    sundials_nvector.h
    sundials_nvector.hpp
    sundials_profiler.h
    sundials_profiler.hpp
    sundials_stepper.h
    sundials_stepper.hpp
    sundials_types_deprecated.h
    sundials_types.h
    sundials_version.h)

if(SUNDIALS_ENABLE_MPI)
  list(APPEND sundials_HEADERS sundials_mpi_errors.h)
endif()

if(SUNDIALS_ENABLE_CUDA)
  list(APPEND sundials_HEADERS sundials_cuda_policies.hpp)
endif()

if(SUNDIALS_ENABLE_HIP)
  list(APPEND sundials_HEADERS sundials_hip_policies.hpp)
endif()

if(SUNDIALS_ENABLE_SYCL)
  list(APPEND sundials_HEADERS sundials_sycl_policies.hpp)
endif()

# If enabled, add the XBraid interface header
if(SUNDIALS_ENABLE_XBRAID)
  list(APPEND sundials_HEADERS sundials_xbraid.h)
endif()

# Add prefix with complete path to the header files
add_prefix(${SUNDIALS_SOURCE_DIR}/include/sundials/ sundials_HEADERS)

set(sundials_SOURCES
    sundatanode/sundatanode_inmem.c
    sundials_adaptcontroller.c
    sundials_adjointcheckpointscheme.c
    sundials_adjointstepper.c
    sundials_band.c
    sundials_cli.c
    sundials_context.c
    sundials_dense.c
    sundials_datanode.c
    sundials_direct.c
    sundials_errors.c
    sundials_domeigestimator.c
    sundials_futils.c
    sundials_hashmap.c
    sundials_iterative.c
    sundials_linearsolver.c
    sundials_logger.c
    sundials_math.c
    sundials_matrix.c
    sundials_memory.c
    sundials_nonlinearsolver.c
    sundials_nvector_senswrapper.c
    sundials_nvector.c
    sundials_stepper.c
    sundials_profiler.c
    sundials_version.c)

if(SUNDIALS_ENABLE_MPI)
  list(APPEND sundials_SOURCES sundials_mpi_errors.c)
endif()

# Add prefix with complete path to the source files
add_prefix(${SUNDIALS_SOURCE_DIR}/src/sundials/ sundials_SOURCES)

if(SUNDIALS_ENABLE_MPI)
  set(_link_mpi_if_needed PUBLIC MPI::MPI_C
                          $<$<LINK_LANGUAGE:CXX>:MPI::MPI_CXX>)
endif()

if(SUNDIALS_ENABLE_PROFILING)
  if(SUNDIALS_ENABLE_CALIPER)
    set(_link_caliper_if_needed PUBLIC caliper)
  endif()
  if(SUNDIALS_ENABLE_ADIAK)
    set(_link_adiak_if_needed PUBLIC adiak::adiak ${CMAKE_DL_LIBS})
  endif()
endif()

# Create a library out of the generic sundials modules
sundials_add_library(
  sundials_core
  SOURCES ${sundials_SOURCES}
  HEADERS ${sundials_HEADERS}
  INCLUDE_SUBDIR sundials
  LINK_LIBRARIES ${_link_mpi_if_needed}
  OUTPUT_NAME sundials_core
  VERSION ${sundialslib_VERSION}
  SOVERSION ${sundialslib_SOVERSION})

# Install private headers
install(
  FILES ${SUNDIALS_SOURCE_DIR}/include/sundials/priv/sundials_context_impl.h
        ${SUNDIALS_SOURCE_DIR}/include/sundials/priv/sundials_errors_impl.h
        ${SUNDIALS_SOURCE_DIR}/include/sundials/priv/sundials_logger_macros.h
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/sundials/priv")

if(SUNDIALS_ENABLE_MPI)
  install(
    FILES
      ${SUNDIALS_SOURCE_DIR}/include/sundials/priv/sundials_mpi_errors_impl.h
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/sundials/priv")
endif()

# SUNDIALS_EXPORT macro
include(GenerateExportHeader)
generate_export_header(
  sundials_core BASE_NAME SUNDIALS EXPORT_FILE_NAME
  "${PROJECT_BINARY_DIR}/include/sundials/sundials_export.h")

# Add F2003 module if the interface is enabled
if(SUNDIALS_ENABLE_FORTRAN)
  add_subdirectory("fmod_int${SUNDIALS_INDEX_SIZE}")
endif()
