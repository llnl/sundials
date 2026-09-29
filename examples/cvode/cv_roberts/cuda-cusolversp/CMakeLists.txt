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
set_source_files_properties(cv_roberts_block_cusolversp_batchqr.cu
                            PROPERTIES CUDA_SOURCE_PROPERTY_FORMAT OBJ)
sundials_add_executable(cv_roberts_block_cusolversp_batchqr
                        cv_roberts_block_cusolversp_batchqr.cu)
set_target_properties(cv_roberts_block_cusolversp_batchqr
                      PROPERTIES FOLDER "Examples")
target_link_libraries(
  cv_roberts_block_cusolversp_batchqr
  PRIVATE sundials_sunlinsolcusolversp sundials_cvode sundials_nveccuda
          ${EXE_EXTRA_LINK_LIBS})
if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
  target_link_libraries(cv_roberts_block_cusolversp_batchqr
                        PRIVATE sundials_cvode_fused_cuda)
endif()

sundials_add_test(
  cv_roberts_block_cusolversp_batchqr cv_roberts_block_cusolversp_batchqr
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE cv_roberts_block_cusolversp_batchqr.out
  EXAMPLE_TYPE develop)

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(
    FILES cv_roberts_block_cusolversp_batchqr.cu
          cv_roberts_block_cusolversp_batchqr.out
    DESTINATION
      ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cv_roberts/cuda-cusolversp)
endif()
