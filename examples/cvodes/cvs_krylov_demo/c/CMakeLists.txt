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
  cvsKrylovDemo_ls cvsKrylovDemo_ls.c
  INSTALL_FILES cvsKrylovDemo_ls_0_1.out cvsKrylovDemo_ls_1.out
                cvsKrylovDemo_ls_2.out)
target_link_libraries(cvsKrylovDemo_ls PRIVATE ${CVODES_C_LIBS}
                                              sundials_nvecserial
                                              ${EXE_EXTRA_LINK_LIBS})
foreach(test_args 1 2)
  sundials_add_test(
    cvsKrylovDemo_ls_${test_args} cvsKrylovDemo_ls
    TEST_ARGS ${test_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE cvsKrylovDemo_ls_${test_args}.out
    EXAMPLE_TYPE develop)
endforeach()
if(SUNDIALS_ENABLE_MONITORING)
  sundials_add_test(
    cvsKrylovDemo_ls_0_1 cvsKrylovDemo_ls
    TEST_ARGS 0 1
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE cvsKrylovDemo_ls_0_1.out
    EXAMPLE_TYPE develop)
endif()

sundials_add_example(cvsKrylovDemo_prec cvsKrylovDemo_prec.c)
target_link_libraries(cvsKrylovDemo_prec PRIVATE ${CVODES_C_LIBS}
                                                sundials_nvecserial
                                                ${EXE_EXTRA_LINK_LIBS})
