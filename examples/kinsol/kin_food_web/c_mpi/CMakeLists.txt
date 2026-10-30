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
foreach(example kin_food_web_kry_bbd_p kin_food_web_kry_p)
  sundials_add_example(
    ${example} SOURCES ${example}.c
    MPI_NPROCS 4
    EXAMPLE_TYPE exclude-single
    LINK_LIBRARIES PRIVATE MPI::MPI_C ${KINSOL_C_LIBS} sundials_nvecparallel
                   ${EXE_EXTRA_LINK_LIBS})
endforeach()
