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
# examples/idas level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# C examples
if(SUNDIALS_ENABLE_C_EXAMPLES)
  add_subdirectory(serial)
  if(SUNDIALS_ENABLE_OPENMP)
    add_subdirectory(C_openmp)
  endif()
  if(SUNDIALS_ENABLE_MPI)
    add_subdirectory(parallel)
  endif()
endif()

if(SUNDIALS_ENABLE_FORTRAN_EXAMPLES)
  add_subdirectory(F2003_serial)
endif()
