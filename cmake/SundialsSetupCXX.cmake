# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ UMBC
#                Cody J. Balos @ LLNL
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
# C++-related tests for SUNDIALS CMake-based configuration.
# ---------------------------------------------------------------

enable_language(CXX)
set(CXX_FOUND TRUE)

# ---------------------------------------------------------------
# Option to specify the C++ standard SUNDIALS will use. Defined here so it is
# set in the same configuration pass as the C++ compiler and related options.
# ---------------------------------------------------------------

# Do not allow decaying to previous standards -- generates error if the standard
# is not supported
sundials_option(CMAKE_CXX_STANDARD_REQUIRED BOOL "Require C++ standard version"
                ON)

if(SUNDIALS_ENABLE_PYTHON
   OR ENABLE_SYCL
   OR ENABLE_GINKGO)
  set(DOCSTR "The C++ standard to use if C++ is enabled (17, 20, 23)")
  sundials_option(CMAKE_CXX_STANDARD STRING "${DOCSTR}" "17" OPTIONS "17;20;23")
else()
  set(DOCSTR "The C++ standard to use if C++ is enabled (14, 17, 20, 23)")
  sundials_option(CMAKE_CXX_STANDARD STRING "${DOCSTR}" "14"
                  OPTIONS "14;17;20;23")
endif()
message(STATUS "CXX standard set to ${CMAKE_CXX_STANDARD}")

set(DOCSTR "Enable C++ compiler specific extensions")
sundials_option(CMAKE_CXX_EXTENSIONS BOOL "${DOCSTR}" ON)
message(STATUS "C++ extensions set to ${CMAKE_CXX_EXTENSIONS}")

# Python interface code requires C++17
if(SUNDIALS_ENABLE_PYTHON AND (CMAKE_CXX_STANDARD LESS "17"))
  message(
    SEND_ERROR
      "CMAKE_CXX_STANDARD must be >= 17 because SUNDIALS_ENABLE_PYTHON=ON")
endif()

# SYCL requires C++17
if(ENABLE_SYCL AND (CMAKE_CXX_STANDARD LESS "17"))
  message(FATAL_ERROR "CMAKE_CXX_STANDARD must be >= 17 because ENABLE_SYCL=ON")
endif()

# Ginkgo requires C++17
if(ENABLE_GINKGO AND (CMAKE_CXX_STANDARD LESS "17"))
  message(
    FATAL_ERROR "CMAKE_CXX_STANDARD must be >= 17 because ENABLE_GINKGO=ON")
endif()
