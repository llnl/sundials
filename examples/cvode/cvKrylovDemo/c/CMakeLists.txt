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
sundials_add_executable(cvKrylovDemo_ls cvKrylovDemo_ls.c)
set_target_properties(cvKrylovDemo_ls PROPERTIES FOLDER "Examples")
target_link_libraries(cvKrylovDemo_ls PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})

sundials_add_test(cvKrylovDemo_ls cvKrylovDemo_ls
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvKrylovDemo_ls.out EXAMPLE_TYPE develop)
sundials_add_test(cvKrylovDemo_ls_1 cvKrylovDemo_ls TEST_ARGS 1
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvKrylovDemo_ls_1.out EXAMPLE_TYPE develop)
sundials_add_test(cvKrylovDemo_ls_2 cvKrylovDemo_ls TEST_ARGS 2
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvKrylovDemo_ls_2.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_MONITORING)
  sundials_add_test(cvKrylovDemo_ls_0_1 cvKrylovDemo_ls TEST_ARGS 0 1
                    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                    ANSWER_FILE cvKrylovDemo_ls_0_1.out EXAMPLE_TYPE develop)
endif()
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvKrylovDemo_ls.c cvKrylovDemo_ls.out cvKrylovDemo_ls_0_1.out
                cvKrylovDemo_ls_1.out cvKrylovDemo_ls_2.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cvKrylovDemo/c)
endif()

sundials_add_executable(cvKrylovDemo_prec cvKrylovDemo_prec.c)
set_target_properties(cvKrylovDemo_prec PROPERTIES FOLDER "Examples")
target_link_libraries(cvKrylovDemo_prec PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvKrylovDemo_prec cvKrylovDemo_prec ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvKrylovDemo_prec.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvKrylovDemo_prec.c cvKrylovDemo_prec.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cvKrylovDemo/c)
endif()
