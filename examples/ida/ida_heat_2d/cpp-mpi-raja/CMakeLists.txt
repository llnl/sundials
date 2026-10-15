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
if(SUNDIALS_RAJA_BACKENDS MATCHES "CUDA")
  set(_lang CUDA)
else()
  set(_lang CXX)
endif()
if((RAJA_BACKENDS MATCHES "TARGET_OPENMP") OR (RAJA_BACKENDS MATCHES "OPENMP"))
  set(_other_libs OpenMP::OpenMP_CXX)
endif()
set_source_files_properties(idaHeat2D_kry_p_mpiraja.cpp PROPERTIES LANGUAGE
                                                                    ${_lang})
sundials_add_example(
  idaHeat2D_kry_p_mpiraja idaHeat2D_kry_p_mpiraja.cpp
  MPI_NPROCS 4
  EXAMPLE_TYPE develop
  INCLUDE_DIRECTORIES PRIVATE ${MPI_CXX_INCLUDE_DIRS}
  LINK_LIBRARIES PRIVATE ${IDA_C_LIBS} sundials_nvecmpiplusx
                 sundials_nvecraja ${MPI_CXX_LIBRARIES}
                 ${EXE_EXTRA_LINK_LIBS} ${_other_libs})
