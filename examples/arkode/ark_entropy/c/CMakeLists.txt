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
sundials_add_example(ark_conserved_exp_entropy_ark ark_conserved_exp_entropy_ark.c
                      TEST_NAME ark_conserved_exp_entropy_ark_1_0
                      TEST_ARGS 1 0
                      ANSWER_FILE ark_conserved_exp_entropy_ark_1_0.out
                      EXAMPLE_TYPE develop
                      INSTALL_FILES ark_conserved_exp_entropy_ark_1_1.out
                      LINK_LIBRARIES PRIVATE sundials_arkode sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(ark_conserved_exp_entropy_ark_1_1
                  ark_conserved_exp_entropy_ark TEST_ARGS 1 1
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE ark_conserved_exp_entropy_ark_1_1.out
                  EXAMPLE_TYPE develop)

sundials_add_example(ark_conserved_exp_entropy_erk ark_conserved_exp_entropy_erk.c
                      TEST_NAME ark_conserved_exp_entropy_erk_1 TEST_ARGS 1
                      ANSWER_FILE ark_conserved_exp_entropy_erk_1.out
                      EXAMPLE_TYPE develop
                      LINK_LIBRARIES PRIVATE sundials_arkode sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})

sundials_add_example(ark_dissipated_exp_entropy ark_dissipated_exp_entropy.c
                      TEST_NAME ark_dissipated_exp_entropy_1_0 TEST_ARGS 1 0
                      ANSWER_FILE ark_dissipated_exp_entropy_1_0.out
                      EXAMPLE_TYPE develop
                      INSTALL_FILES ark_dissipated_exp_entropy_1_1.out
                      LINK_LIBRARIES PRIVATE sundials_arkode sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(ark_dissipated_exp_entropy_1_1
                  ark_dissipated_exp_entropy TEST_ARGS 1 1
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE ark_dissipated_exp_entropy_1_1.out
                  EXAMPLE_TYPE develop)
