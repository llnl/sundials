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
# CVODE RAJA C++ example
set(CVODE_examples "cv_adv_diff_kry_raja.cpp\;develop")
if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
  set(_fused_stubs sundials_cvode_fused_stubs)
endif()
if(SUNDIALS_RAJA_BACKENDS MATCHES "CUDA")
  set(_lang CUDA)
else()
  set(_lang CXX)
endif()
if((RAJA_BACKENDS MATCHES "TARGET_OPENMP") OR (RAJA_BACKENDS MATCHES "OPENMP"))
  set(_other_libs OpenMP::OpenMP_CXX)
endif()
set_source_files_properties(cv_adv_diff_kry_raja.cpp PROPERTIES LANGUAGE ${_lang})
sundials_add_executable(cv_adv_diff_kry_raja cv_adv_diff_kry_raja.cpp)
set_target_properties(cv_adv_diff_kry_raja PROPERTIES FOLDER "Examples")
target_include_directories(cv_adv_diff_kry_raja PRIVATE .)
target_link_libraries(cv_adv_diff_kry_raja
  PRIVATE sundials_cvode sundials_nvecraja ${_fused_stubs} ${EXE_EXTRA_LINK_LIBS}
          ${_other_libs})
sundials_add_test(cv_adv_diff_kry_raja cv_adv_diff_kry_raja
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cv_adv_diff_kry_raja.out EXAMPLE_TYPE develop)

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cv_adv_diff_kry_raja.cpp cv_adv_diff_kry_raja.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cvAdvDiff/cpp-raja)
endif()
