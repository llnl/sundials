# -----------------------------------------------------------------
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
# -----------------------------------------------------------------
sundials_add_example(
  test_fnvector_fortran_mod
  SOURCES test_fnvector_fortran_mod.f90 fnvector_fortran_mod.f90
  NO_DIFF
  EXAMPLE_TYPE develop
  PROPERTIES
    Fortran_MODULE_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/test_fnvector_fortran_mod.dir
  LINK_LIBRARIES PRIVATE sundials_fcore_mod ${EXE_EXTRA_LINK_LIBS})
