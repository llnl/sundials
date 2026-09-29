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
sundials_add_executable(cv_adv_diff_diag_hip cv_adv_diff_diag_hip.cpp)
set_target_properties(cv_adv_diff_diag_hip PROPERTIES FOLDER "Examples")
target_link_libraries(
  cv_adv_diff_diag_hip PRIVATE hip::device sundials_cvode sundials_nvechip
                               ${EXE_EXTRA_LINK_LIBS})
if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
  target_link_libraries(cv_adv_diff_diag_hip PRIVATE sundials_cvode_fused_hip)
endif()
foreach(example_args "0 0" "0 1" "1 1")
  string(REGEX REPLACE " " "_" test_name cv_adv_diff_diag_hip_${example_args})
  sundials_add_test(
    ${test_name} cv_adv_diff_diag_hip
    TEST_ARGS ${example_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE develop)
endforeach()
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  file(GLOB example_out cv_adv_diff_diag_hip*.out)
  install(FILES cv_adv_diff_diag_hip.cpp ${example_out}
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cv_adv_diff/hip)
endif()

sundials_add_executable(cv_adv_diff_kry_hip cv_adv_diff_kry_hip.cpp)
set_target_properties(cv_adv_diff_kry_hip PROPERTIES FOLDER "Examples")
target_link_libraries(
  cv_adv_diff_kry_hip PRIVATE hip::device sundials_cvode sundials_nvechip
                              ${EXE_EXTRA_LINK_LIBS})
if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
  target_link_libraries(cv_adv_diff_kry_hip PRIVATE sundials_cvode_fused_hip)
endif()
sundials_add_test(
  cv_adv_diff_kry_hip cv_adv_diff_kry_hip
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE cv_adv_diff_kry_hip.out
  EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cv_adv_diff_kry_hip.cpp cv_adv_diff_kry_hip.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cv_adv_diff/hip)
endif()
