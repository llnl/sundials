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
# CVODE MAGMA C++ example
set(CVODE_examples "cv_bruss_batched_magma.cpp\;\;develop")
if(SUNDIALS_MAGMA_BACKENDS MATCHES "CUDA")
  set(vector nveccuda)
  set(cuda_or_hip CUDA)
elseif(SUNDIALS_MAGMA_BACKENDS MATCHES "HIP")
  set(vector nvechip)
  set(cuda_or_hip HIP)
endif()
set_source_files_properties(cv_bruss_batched_magma.cpp
                            PROPERTIES LANGUAGE ${cuda_or_hip})
sundials_add_executable(cv_bruss_batched_magma cv_bruss_batched_magma.cpp)
set_target_properties(cv_bruss_batched_magma PROPERTIES FOLDER "Examples")
target_link_libraries(cv_bruss_batched_magma
  PRIVATE sundials_cvode sundials_${vector} sundials_sunlinsolmagmadense
          ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cv_bruss_batched_magma cv_bruss_batched_magma
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cv_bruss_batched_magma.out EXAMPLE_TYPE develop)

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cv_bruss_batched_magma.cpp cv_bruss_batched_magma.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cv_brusselator/cpp-magma)
endif()
