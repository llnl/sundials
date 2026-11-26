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
# Setup the HIP language and libraries.
# ---------------------------------------------------------------

if(NOT DEFINED ROCM_PATH)
  if(NOT DEFINED ENV{ROCM_PATH})
    set(ROCM_PATH
        "/opt/rocm/"
        CACHE PATH "Path to which ROCm has been installed")
  else()
    set(ROCM_PATH
        "$ENV{ROCM_PATH}"
        CACHE PATH "Path to which ROCm has been installed")
  endif()
endif()

if(NOT DEFINED HIP_PATH)
  if(NOT DEFINED ENV{HIP_PATH})
    set(HIP_PATH
        "/opt/rocm/hip"
        CACHE PATH "Path to which HIP has been installed")
  else()
    set(HIP_PATH
        "$ENV{HIP_PATH}"
        CACHE PATH "Path to which HIP has been installed")
  endif()
endif()

if(NOT DEFINED HIP_PLATFORM)
  if(NOT DEFINED ENV{HIP_PLATFORM})
    set(HIP_PLATFORM
        "amd"
        CACHE STRING "HIP platform (amd, nvidia)")
  else()
    set(HIP_PLATFORM
        "$ENV{HIP_PLATFORM}"
        CACHE STRING "HIP platform (amd, nvidia)")
  endif()
endif()

# ===============================================================
# Enable the HIP language
# ===============================================================

# A C++ compiler that is not Clang based cannot compile HIP sources, so they are
# compiled with the HIP language instead. SUNDIALS_HIP_LANGUAGE is the language
# of HIP sources and SUNDIALS_HIP_TARGET the HIP target to link: hip::device
# compiles every C++ source as HIP, so it requires a HIP compiler as C++
# compiler.
if(HIP_PLATFORM STREQUAL "amd" AND NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  if(CMAKE_VERSION VERSION_LESS 3.21)
    message(FATAL_ERROR "The HIP language requires CMake 3.21 or newer")
  endif()

  if(NOT DEFINED CMAKE_HIP_ARCHITECTURES)
    if(AMDGPU_TARGETS)
      set(CMAKE_HIP_ARCHITECTURES "${AMDGPU_TARGETS}")
    elseif(GPU_TARGETS)
      set(CMAKE_HIP_ARCHITECTURES "${GPU_TARGETS}")
    endif()
  endif()

  sundials_option(CMAKE_HIP_STANDARD_REQUIRED BOOL
                  "Require HIP standard version" ON)

  set(DOCSTR "The HIP standard to use if HIP is enabled (14, 17, 20, 23)")
  sundials_option(CMAKE_HIP_STANDARD STRING "${DOCSTR}" "${CMAKE_CXX_STANDARD}"
                  OPTIONS "14;17;20;23")
  message(STATUS "HIP standard set to ${CMAKE_HIP_STANDARD}")

  enable_language(HIP)
  message(STATUS "HIP language compiler: ${CMAKE_HIP_COMPILER}")
  message(STATUS "HIP architectures: ${CMAKE_HIP_ARCHITECTURES}")

  set(SUNDIALS_HIP_LANGUAGE HIP)
  set(SUNDIALS_HIP_TARGET hip::host)
else()
  set(SUNDIALS_HIP_LANGUAGE CXX)
  set(SUNDIALS_HIP_TARGET hip::device)
endif()

# ===============================================================
# Find HIP
# ===============================================================

# Set CMAKE_PREFIX_PATH as the hip-config.cmake has some find_package calls
# which don't have the proper path set (not sure if this is a bug or
# intentional), so without this they will fail even if we provide the PATH
# option to find_package(HIP).
set(CMAKE_PREFIX_PATH "${ROCM_PATH};${HIP_PATH}")
find_package(HIP REQUIRED)

if("${HIP_COMPILER}" STREQUAL "hcc")
  message(FATAL_ERROR "Deprecated HCC compiler is not supported"
                      "Please update ROCm")
endif()

message(STATUS "HIP version:      ${HIP_VERSION}")
message(STATUS "HIP platform:     ${HIP_PLATFORM}")
message(STATUS "HIP compiler:     ${HIP_COMPILER}")
message(STATUS "HIP linker:       ${CMAKE_CXX_LINK_EXECUTABLE}")
message(STATUS "AMD targets:      ${AMDGPU_TARGETS}")
