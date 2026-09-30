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
  cv_roberts_block_klu cv_roberts_block_klu.c ADD_TEST
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE cv_roberts_block_klu.out
  EXAMPLE_TYPE develop
  INSTALL_FILES cv_roberts_block_klu.out)
set_target_properties(cv_roberts_block_klu PROPERTIES FOLDER "Examples")
target_link_libraries(
  cv_roberts_block_klu PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                               sundials_sunlinsolklu ${EXE_EXTRA_LINK_LIBS})

sundials_add_example(
  cv_roberts_klu cv_roberts_klu.c ADD_TEST
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE cv_roberts_klu.out
  EXAMPLE_TYPE develop
  INSTALL_FILES cv_roberts_klu.out)
set_target_properties(cv_roberts_klu PROPERTIES FOLDER "Examples")
target_link_libraries(
  cv_roberts_klu PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                         sundials_sunlinsolklu ${EXE_EXTRA_LINK_LIBS})
