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
# CVODE Kokkos batched Brusselator example
foreach(backend ${KOKKOS_EXAMPLES_BACKENDS})
  set(example_target cv_bruss_batched_kokkos.${backend})
  if(NOT TARGET ${example_target})
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${example_target}.out")
      message(STATUS "Skipping ${example_target}: no answer file is available")
      continue()
    endif()
    sundials_add_example(
      ${example_target} cv_bruss_batched_kokkos.cpp
      EXTRA_ARGS --kokkos-disable-warnings
      INSTALL_FILES
        "${PROJECT_SOURCE_DIR}/examples/utilities/example_utilities.hpp"
      LINK_LIBRARIES
        PRIVATE
        ${CVODE_C_LIBS}
        sundials_nveckokkos
        sundials_sunmatrixkokkosdense
        sundials_sunlinsolkokkosdense
        ${EXE_EXTRA_LINK_LIBS}
      INCLUDE_DIRECTORIES PRIVATE ${PROJECT_SOURCE_DIR}/examples/utilities
      COMPILE_DEFINITIONS PRIVATE USE_${backend})
  endif()
endforeach()
