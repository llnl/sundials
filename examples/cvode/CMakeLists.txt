# ---------------------------------------------------------------
# Programmer(s): David J. Gardner, Slaven Peles, and
#                Shelby Lockhart @ LLNL
# ---------------------------------------------------------------
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
# ---------------------------------------------------------------
# examples/cvode level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# Keep the common CVODE library list available to the C variants below.
set(CVODE_C_LIBS sundials_cvode)
if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
  list(APPEND CVODE_C_LIBS sundials_cvode_fused_stubs)
endif()

# Examples are grouped by problem. Each problem directory contains variants
# grouped by implementation language and optional TPLs.
add_subdirectory(cvAdvDiff)
add_subdirectory(cvAnalytic)
add_subdirectory(cvAnalyticSys)
add_subdirectory(cvBrusselator)
add_subdirectory(cvDiag)
add_subdirectory(cvDirectDemo)
add_subdirectory(cvDisc)
add_subdirectory(cvDiurnal)
add_subdirectory(cvHeat2D)
add_subdirectory(cvKPR)
add_subdirectory(cvKrylovDemo)
add_subdirectory(cvParticle)
add_subdirectory(cvPendulum)
add_subdirectory(cvPetsc)
add_subdirectory(cvRoberts)
add_subdirectory(cvRocket)
add_subdirectory(cvVdp)

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  # Install the reorganized source CMake files below the generated CVODE root
  # CMakeLists.txt. The installed root provides the small wrappers needed by
  # these in-tree files and selects variants using the configured SUNDIALS
  # options.
  file(GLOB_RECURSE _cvode_cmakelists RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
       "cv*/CMakeLists.txt")
  foreach(_cmakelist ${_cvode_cmakelists})
    get_filename_component(_cmakelist_dir ${_cmakelist} DIRECTORY)
    install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/${_cmakelist}
            DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/${_cmakelist_dir})
  endforeach()

  set(CVODE_LANGUAGES C)
  if(SUNDIALS_ENABLE_CXX_EXAMPLES)
    list(APPEND CVODE_LANGUAGES CXX)
  endif()
  if(SUNDIALS_ENABLE_FORTRAN_EXAMPLES)
    list(APPEND CVODE_LANGUAGES Fortran)
  endif()
  if(SUNDIALS_ENABLE_CUDA_EXAMPLES AND SUNDIALS_ENABLE_CUDA)
    list(APPEND CVODE_LANGUAGES CUDA)
  endif()
  list2string(CVODE_LANGUAGES CVODE_LANGUAGES_STRING)

  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_cvode_ex.in
    ${PROJECT_BINARY_DIR}/examples/cvode/CMakeLists.txt @ONLY)
  install(FILES ${PROJECT_BINARY_DIR}/examples/cvode/CMakeLists.txt
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode)

  # The installed CVODE tree is one project rooted at cvode/, so use the
  # installation smoke-test helper with the root directory as its test path.
  sundials_add_test_install(cvode .)

  install(FILES README.md DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode)
endif()
