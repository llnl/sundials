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
# CVODE Ginkgo batched Brusselator example
foreach(backend REF OMP CUDA HIP SYCL)
  if(NOT SUNDIALS_GINKGO_BACKENDS MATCHES "${backend}")
    continue()
  endif()

  set(float_precision default)
  if(backend MATCHES "CUDA")
    set(vector nveccuda)
    set(language CUDA)
    set(float_precision 4)
  elseif(backend MATCHES "HIP")
    set(vector nvechip)
    set(language CXX)
  elseif(backend MATCHES "SYCL")
    set(vector nvecsycl)
    set(language CXX)
  elseif(backend MATCHES "OMP")
    set(vector nvecopenmp)
    set(language CXX)
  else()
    set(vector nvecserial)
    set(language CXX)
  endif()

  set_source_files_properties(cv_bruss_batched_ginkgo.cpp
                              PROPERTIES LANGUAGE ${language})
  set(example_target cv_bruss_batched_ginkgo.${backend})
  set(example_options)
  if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${example_target}.out")
    message(
      STATUS
        "No answer file for ${example_target}; running without output comparison"
    )
    list(APPEND example_options NODIFF)
  endif()
  sundials_add_example(
    ${example_target} cv_bruss_batched_ginkgo.cpp ${example_options}
    FLOAT_PRECISION ${float_precision}
    INSTALL_FILES
      cv_bruss_batched_ginkgo.hpp
      "${PROJECT_SOURCE_DIR}/examples/utilities/example_utilities.hpp"
    LINK_LIBRARIES PRIVATE ${CVODE_C_LIBS} sundials_${vector} Ginkgo::ginkgo
                   ${EXE_EXTRA_LINK_LIBS}
    INCLUDE_DIRECTORIES PRIVATE ${PROJECT_SOURCE_DIR}/examples/utilities
    COMPILE_DEFINITIONS PRIVATE USE_${backend})
endforeach()
