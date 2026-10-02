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
  include(SundialsInstallExampleCMake)

  # Install each variant as an independent CMake project.
  get_property(
    _cvode_problem_directories
    DIRECTORY
    PROPERTY SUBDIRECTORIES)
  set(_cvode_example_directories)
  foreach(_problem_dir ${_cvode_problem_directories})
    get_property(
      _variant_dirs
      DIRECTORY ${_problem_dir}
      PROPERTY SUBDIRECTORIES)
    list(APPEND _cvode_example_directories ${_variant_dirs})
  endforeach()

  foreach(_example_dir ${_cvode_example_directories})
    file(RELATIVE_PATH _example_rel_dir ${CMAKE_CURRENT_SOURCE_DIR}
         ${_example_dir})
    sundials_install_example_cmake(
      "${_example_dir}"
      "${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/${_example_rel_dir}")
  endforeach()

  sundials_add_test_install(cvode cv_adv_diff/c EXECUTABLE cv_adv_diff_bnd)
endif()
