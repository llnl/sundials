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
if(SUNDIALS_MAGMA_BACKENDS MATCHES "CUDA")
  set(vector nveccuda)
  set(language CUDA)
elseif(SUNDIALS_MAGMA_BACKENDS MATCHES "HIP")
  set(vector nvechip)
  set(language CXX)
endif()
set_source_files_properties(cv_bruss_batched_magma.cpp PROPERTIES LANGUAGE
                                                                  ${language})
sundials_add_example(cv_bruss_batched_magma cv_bruss_batched_magma.cpp)
target_link_libraries(
  cv_bruss_batched_magma
  PRIVATE sundials_cvode sundials_${vector} sundials_sunlinsolmagmadense
          ${EXE_EXTRA_LINK_LIBS})
