# ---------------------------------------------------------------
# Programmer(s): Daniel R. Reynolds @ UMBC
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
# unit_tests/arkode level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# C unit tests
add_subdirectory(C_serial)

# C++ unit tests
if(CXX_FOUND)
  add_subdirectory(CXX_serial)
  if(SUNDIALS_TEST_ENABLE_GTEST)
    add_subdirectory(gtest)
  endif()
  if(ENABLE_MPI AND MPI_CXX_FOUND)
    add_subdirectory(CXX_parallel)
  endif()
endif()

# F2003 unit tests
if(BUILD_FORTRAN_MODULE_INTERFACE)
  add_subdirectory(F2003_serial)
endif()
