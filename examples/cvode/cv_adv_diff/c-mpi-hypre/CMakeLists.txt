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
  cv_adv_diff_non_ph cv_adv_diff_non_ph.c ADD_TEST
  MPI_NPROCS 2
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE cv_adv_diff_non_ph.out
  EXAMPLE_TYPE develop
  INSTALL_FILES cv_adv_diff_non_ph.c cv_adv_diff_non_ph.out)
set_target_properties(cv_adv_diff_non_ph PROPERTIES FOLDER "Examples")
target_link_libraries(
  cv_adv_diff_non_ph PRIVATE ${CVODE_C_LIBS} sundials_nvecparhyp
                             ${EXE_EXTRA_LINK_LIBS})
