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
foreach(_example ark_diag_kry_bbd_f2003 ark_diag_non_f2003)
  sundials_add_example(${_example} ${_example}.f90 MPI_NPROCS 4
                        FLOAT_PRECISION 2 EXAMPLE_TYPE develop
                        LINK_LIBRARIES PRIVATE MPI::MPI_Fortran
                        sundials_arkode sundials_farkode_mod
                        sundials_nvecparallel sundials_fnvecparallel_mod
                        sundials_nvecmpiplusx sundials_fnvecmpiplusx_mod
                        ${EXE_EXTRA_LINK_LIBS})
endforeach()
