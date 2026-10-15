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
  idasKrylovDemo_ls idasKrylovDemo_ls.c
  INSTALL_FILES idasKrylovDemo_ls_1.out idasKrylovDemo_ls_2.out)
target_link_libraries(idasKrylovDemo_ls PRIVATE ${IDAS_C_LIBS}
                                                sundials_nvecserial
                                                ${EXE_EXTRA_LINK_LIBS})
foreach(test_args 1 2)
  sundials_add_test(
    idasKrylovDemo_ls_${test_args} idasKrylovDemo_ls
    TEST_ARGS ${test_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE idasKrylovDemo_ls_${test_args}.out
    EXAMPLE_TYPE develop)
endforeach()
