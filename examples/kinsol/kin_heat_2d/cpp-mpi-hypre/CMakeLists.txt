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
  kin_heat2D_nonlin_hypre_pfmg kin_heat2D_nonlin_hypre_pfmg.cpp
  TEST_NAME kin_heat2D_nonlin_hypre_pfmg_--np_2_2
  TEST_ARGS --np 2 2
  MPI_NPROCS 4
  INSTALL_FILES kin_heat2D_nonlin_hypre_pfmg.hpp
  LINK_LIBRARIES PRIVATE MPI::MPI_CXX ${KINSOL_C_LIBS} sundials_nvecparallel
                 SUNDIALS::HYPRE ${EXE_EXTRA_LINK_LIBS})
