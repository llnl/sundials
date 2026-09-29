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
add_subdirectory(cv_adv_diff)
add_subdirectory(cv_analytic)
add_subdirectory(cv_analytic_sys)
add_subdirectory(cv_brusselator)
add_subdirectory(cv_diag)
add_subdirectory(cv_direct_demo)
add_subdirectory(cv_disc)
add_subdirectory(cv_diurnal)
add_subdirectory(cv_heat_2d)
add_subdirectory(cv_kpr)
add_subdirectory(cv_krylov_demo)
add_subdirectory(cv_particle)
add_subdirectory(cv_pendulum)
add_subdirectory(cv_petsc)
add_subdirectory(cv_roberts)
add_subdirectory(cv_rocket)
add_subdirectory(cv_vdp)

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  # Install each variant as an independent CMake project. The configured
  # CMakeLists.txt supplies the installed-package setup and then includes the
  # variant-specific build rules.
  file(
    GLOB_RECURSE _cvode_cmakelists
    RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    "cv*/*/CMakeLists.txt")

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
    ${PROJECT_BINARY_DIR}/examples/cvode/CMakeLists.installed.txt @ONLY)

  foreach(_cmakelist ${_cvode_cmakelists})
    get_filename_component(_cmakelist_dir ${_cmakelist} DIRECTORY)
    install(
      FILES ${CMAKE_CURRENT_SOURCE_DIR}/${_cmakelist}
      DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/${_cmakelist_dir}
      RENAME cvode-example.cmake)
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/cvode/CMakeLists.installed.txt
      DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/${_cmakelist_dir}
      RENAME CMakeLists.txt)
  endforeach()

  sundials_add_test_install(cvode cv_adv_diff/c EXECUTABLE cv_adv_diff_bnd)
endif()
