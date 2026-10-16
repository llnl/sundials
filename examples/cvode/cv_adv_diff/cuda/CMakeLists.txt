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
set_source_files_properties(cv_adv_diff_diag_cuda.cu
                            PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
sundials_add_example(
  cv_adv_diff_diag_cuda cv_adv_diff_diag_cuda.cu NO_TEST
  INSTALL_FILES cv_adv_diff_diag_cuda_0_0.out cv_adv_diff_diag_cuda_0_1.out
                cv_adv_diff_diag_cuda_1_1.out)
target_link_libraries(
  cv_adv_diff_diag_cuda PRIVATE sundials_cvode sundials_nveccuda
                                ${EXE_EXTRA_LINK_LIBS})
if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
  target_link_libraries(cv_adv_diff_diag_cuda PRIVATE sundials_cvode_fused_cuda)
endif()

foreach(example_args "0 0" "0 1" "1 1")
  string(REGEX REPLACE " " "_" test_name cv_adv_diff_diag_cuda_${example_args})
  sundials_add_test(
    ${test_name} cv_adv_diff_diag_cuda
    TEST_ARGS ${example_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${test_name}.out
    EXAMPLE_TYPE develop)
endforeach()

set_source_files_properties(cv_adv_diff_kry_cuda.cu
                            PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
sundials_add_example(cv_adv_diff_kry_cuda cv_adv_diff_kry_cuda.cu)
target_link_libraries(
  cv_adv_diff_kry_cuda PRIVATE sundials_cvode sundials_nveccuda
                               ${EXE_EXTRA_LINK_LIBS})
if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
  target_link_libraries(cv_adv_diff_kry_cuda PRIVATE sundials_cvode_fused_cuda)
endif()
