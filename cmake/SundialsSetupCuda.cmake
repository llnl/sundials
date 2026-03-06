# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# Setup the CUDA language and CUDA libraries.
# ---------------------------------------------------------------

# ===============================================================
# Configure options needed prior to enabling the CUDA language
# ===============================================================

if(NOT CMAKE_CUDA_HOST_COMPILER)
  # If a user did not provide the host compiler, then we assume that they want
  # to use the CXX compiler that was set.
  set(CMAKE_CUDA_HOST_COMPILER
      ${CMAKE_CXX_COMPILER}
      CACHE FILEPATH "NVCC host compiler")
endif()

# ===============================================================
# Configure the CUDA flags
# ===============================================================

# Do not allow decaying to previous standards -- generates error if the standard
# is not supported
sundials_option(CMAKE_CUDA_STANDARD_REQUIRED BOOL
                "Require C++ standard version" ON)

set(DOCSTR "The CUDA standard to use if CUDA is enabled (14, 17, 20)")
sundials_option(CMAKE_CUDA_STANDARD STRING "${DOCSTR}" "${CMAKE_CXX_STANDARD}"
                OPTIONS "14;17;20")
message(STATUS "CUDA standard set to ${CMAKE_CUDA_STANDARD}")

set(CMAKE_CUDA_FLAGS
    "${CMAKE_CUDA_FLAGS} --expt-extended-lambda --expt-relaxed-constexpr")

if((CMAKE_CXX_COMPILER_ID MATCHES GNU)
   OR (CMAKE_CXX_COMPILER_ID MATCHES Clang)
   AND (CMAKE_SYSTEM_PROCESSOR MATCHES ppc64le))
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-mno-float128 _hasflag)
  if(_hasflag)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=-mno-float128")
  endif()
endif()

# ===============================================================
# Prohibit CUDA interface when using extended precision.
# ===============================================================

if(SUNDIALS_PRECISION MATCHES "EXTENDED")
  message(
    FATAL_ERROR "CUDA interfaces are incompatible with extended precision.")
endif()

# ===============================================================
# Enable CUDA lang and find the CUDA libraries.
# ===============================================================

enable_language(CUDA)
set(CUDA_FOUND TRUE)

find_package(CUDAToolkit REQUIRED)

# Show CUDA flags
mark_as_advanced(CLEAR CMAKE_CUDA_FLAGS)

# ===============================================================
# Print out information about CUDA.
# ===============================================================

message(STATUS "CUDA Toolkit Version:       ${CUDAToolkit_VERSION}")
message(STATUS "CUDA Architectures:         ${CMAKE_CUDA_ARCHITECTURES}")
message(STATUS "CUDA Compiler:              ${CMAKE_CUDA_COMPILER}")
message(STATUS "CUDA Host Compiler:         ${CMAKE_CUDA_HOST_COMPILER}")
message(STATUS "CUDA Toolkit Includes:      ${CUDAToolkit_INCLUDE_DIRS}")
message(STATUS "CUDA Library Directory:     ${CUDAToolkit_LIBRARY_DIR}")
message(STATUS "CUDA Compile Flags:         ${CMAKE_CUDA_FLAGS}")
message(STATUS "CUDA Link Flags:            ${CMAKE_CUDA_LINK_FLAGS}")
message(STATUS "CUDA Link Executable:       ${CMAKE_CUDA_LINK_EXECUTABLE}")
message(
  STATUS "CUDA Separable Compilation: ${CMAKE_CUDA_SEPARABLE_COMPILATION}")

# ===============================================================
# Configure compiler for installed examples
# ===============================================================

if(ENABLE_MPI)
  set(_EXAMPLES_CUDA_HOST_COMPILER
      "${MPI_CXX_COMPILER}"
      CACHE INTERNAL "${lang} compiler for installed examples")
else()
  set(_EXAMPLES_CUDA_HOST_COMPILER
      "${CMAKE_CUDA_HOST_COMPILER}"
      CACHE INTERNAL "${lang} compiler for installed examples")
endif()
