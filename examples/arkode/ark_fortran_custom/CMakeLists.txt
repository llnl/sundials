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
if(SUNDIALS_ENABLE_FORTRAN_EXAMPLES AND SUNDIALS_INDEX_SIZE MATCHES "64")
  add_subdirectory(fortran_complex)
  add_subdirectory(fortran_brusselator_1d)
  add_subdirectory(fortran_fnvector_fortran_mod)
  add_subdirectory(fortran_fnvector_complex_mod)
  add_subdirectory(fortran_fsunmatrix_fortran_mod)
  add_subdirectory(fortran_fsunlinsol_fortran_mod)
endif()
