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
if(SUNDIALS_INDEX_SIZE MATCHES "64")
  sundials_add_example(
    ida_heat_2d_kry_bbd_f2003
    SOURCES ida_heat_2d_kry_bbd_f2003.f90
    MPI_NPROCS 4
    EXAMPLE_TYPE develop
    FLOAT_PRECISION 2
    INSTALL_FILES ida_heat_2d_kry_bbd_f2003_4.out
    LINK_LIBRARIES
      PRIVATE
      MPI::MPI_Fortran
      ${IDA_C_LIBS}
      sundials_fida_mod
      sundials_nvecparallel
      sundials_fnvecparallel_mod
      ${EXE_EXTRA_LINK_LIBS}
    PROPERTIES
      Fortran_MODULE_DIRECTORY
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/ida_heat_2d_kry_bbd_f2003.dir)
endif()
