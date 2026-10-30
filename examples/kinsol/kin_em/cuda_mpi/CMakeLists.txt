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
  kin_em_mpicuda
  SOURCES kin_em_mpicuda.cu
  MPI_NPROCS 2
  INSTALL_FILES kin_em_mpicuda.hpp
  INCLUDE_DIRECTORIES PRIVATE ${MPI_CXX_INCLUDE_DIRS}
  LINK_LIBRARIES PRIVATE ${KINSOL_C_LIBS} sundials_nvecmpiplusx
                 sundials_nveccuda ${MPI_CXX_LIBRARIES} ${EXE_EXTRA_LINK_LIBS})
