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
  ark_adr_1d_extsts
  SOURCES ark_adr_1d_extsts.cpp
  EXAMPLE_TYPE develop
  INSTALL_FILES
    ark_adr_1d_extsts.hpp
    ark_adr_1d_extsts.out
    ark_adr1d_extsts_arkode.use_extsts_error_estimate_1.out
  INCLUDE_DIRECTORIES PRIVATE ${PROJECT_SOURCE_DIR}/examples/utilities
  LINK_LIBRARIES PRIVATE sundials_arkode sundials_nvecserial
                 ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(
  ark_adr_1d_extsts ark_adr_1d_extsts
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE ark_adr_1d_extsts.out)
sundials_add_test(
  ark_adr1d_extsts_arkode.use_extsts_error_estimate_1
  ark_adr_1d_extsts
  TEST_ARGS arkode.use_extsts_error_estimate 1
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE
    ark_adr1d_extsts_arkode.use_extsts_error_estimate_1.out
)
