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
  kin_heat2D_nonlin_p kin_heat2D_nonlin_p.cpp
  TEST_NAME kin_heat2D_nonlin_p_--np_2_2
  TEST_ARGS --np 2 2
  MPI_NPROCS 4
  EXAMPLE_TYPE exclude-single
  INSTALL_FILES kin_heat2D_nonlin_p.hpp
  LINK_LIBRARIES PRIVATE ${KINSOL_C_LIBS} sundials_nvecparallel
                 sundials_nvecmpiplusx sundials_nvecserial MPI::MPI_CXX
                 ${EXE_EXTRA_LINK_LIBS})
