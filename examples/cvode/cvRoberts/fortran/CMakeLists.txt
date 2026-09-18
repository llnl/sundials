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
set(example cv_roberts_dns_f2003)
set(CMAKE_Fortran_MODULE_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${example}.dir)
sundials_add_executable(${example} ${example}.f90)
set_target_properties(${example} PROPERTIES FOLDER "Examples")
sundials_add_test(
  ${example} ${example}
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE ${example}.out
  EXAMPLE_TYPE develop)
target_link_libraries(${example} sundials_cvode sundials_fcvode_mod
                      ${EXE_EXTRA_LINK_LIBS})
if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
  target_link_libraries(${example} sundials_cvode_fused_stubs)
endif()
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES ${example}.f90 ${example}.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cvRoberts/fortran)
endif()

set(example cv_roberts_dns_constraints_f2003)
set(CMAKE_Fortran_MODULE_DIRECTORY
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${example}.dir)
sundials_add_executable(${example} ${example}.f90)
set_target_properties(${example} PROPERTIES FOLDER "Examples")
sundials_add_test(
  ${example} ${example}
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE ${example}.out
  EXAMPLE_TYPE develop)
target_link_libraries(${example} sundials_cvode sundials_fcvode_mod
                      ${EXE_EXTRA_LINK_LIBS})
if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
  target_link_libraries(${example} sundials_cvode_fused_stubs)
endif()
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES ${example}.f90 ${example}.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cvRoberts/fortran)
endif()
