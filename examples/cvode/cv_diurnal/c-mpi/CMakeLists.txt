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
  cv_diurnal_kry_bbd_p cv_diurnal_kry_bbd_p.c
  MPI_NPROCS 4
  EXAMPLE_TYPE exclude-single
  LINK_LIBRARIES PRIVATE MPI::MPI_C ${CVODE_C_LIBS} sundials_nvecparallel
                 ${EXE_EXTRA_LINK_LIBS})

sundials_add_example(
  cv_diurnal_kry_p cv_diurnal_kry_p.c
  MPI_NPROCS 4
  EXAMPLE_TYPE exclude-single
  LINK_LIBRARIES PRIVATE MPI::MPI_C ${CVODE_C_LIBS} sundials_nvecparallel
                 ${EXE_EXTRA_LINK_LIBS})
