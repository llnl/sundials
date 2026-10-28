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
sundials_add_example(cv_adv_diff_petsc cv_adv_diff_petsc.c MPI_NPROCS 4)
target_link_libraries(
  cv_adv_diff_petsc
  PRIVATE MPI::MPI_C ${CVODE_C_LIBS} sundials_nvecpetsc
          sundials_sunnonlinsolpetscsnes ${EXE_EXTRA_LINK_LIBS})
