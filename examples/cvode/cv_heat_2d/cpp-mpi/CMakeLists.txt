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
# CVODE parallel C++ heat equation example
sundials_add_example(
  cv_heat_2d_p cv_heat_2d_p.cpp ADD_TEST
  TEST_ARGS --np 2 2
  MPI_NPROCS 4
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE cv_heat_2d_p_--np_2_2.out
  EXAMPLE_TYPE develop
  INSTALL_FILES cv_heat_2d_p.cpp cv_heat_2d_p_--np_2_2.out plot_heat_2d_p.py)
set_target_properties(cv_heat_2d_p PROPERTIES FOLDER "Examples")
target_link_libraries(
  cv_heat_2d_p PRIVATE MPI::MPI_CXX ${CVODE_C_LIBS} sundials_nvecparallel
                       ${EXE_EXTRA_LINK_LIBS})
