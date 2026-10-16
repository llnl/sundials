# ---------------------------------------------------------------
# Programmer(s): David J. Gardner and Slaven Peles @ LLNL
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
# examples/arkode level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

# Examples are grouped by problem. Each problem directory contains variants
# grouped by implementation language and optional third-party libraries.
add_subdirectory(ark_advection_diffusion_reaction)
add_subdirectory(ark_analytic)
add_subdirectory(ark_brusselator)
add_subdirectory(ark_butcher)
add_subdirectory(ark_diffusion_reaction)
add_subdirectory(ark_diurnal)
add_subdirectory(ark_entropy)
add_subdirectory(ark_fortran_custom)
add_subdirectory(ark_harmonic)
add_subdirectory(ark_heat_1d)
add_subdirectory(ark_heat_2d)
add_subdirectory(ark_kepler)
add_subdirectory(ark_kpr)
add_subdirectory(ark_krylov_demo)
add_subdirectory(ark_lotka_volterra)
add_subdirectory(ark_mri_step)
add_subdirectory(ark_pendulum)
add_subdirectory(ark_petsc)
add_subdirectory(ark_robertson)
add_subdirectory(ark_sod)

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  include(SundialsInstallExampleCMake)

  # Install each language/TPL variant as an independent CMake project, as is
  # done for the CVODE examples.
  get_property(
    _arkode_problem_directories
    DIRECTORY
    PROPERTY SUBDIRECTORIES)
  set(_arkode_example_directories)
  foreach(_problem_dir ${_arkode_problem_directories})
    get_property(
      _variant_dirs
      DIRECTORY ${_problem_dir}
      PROPERTY SUBDIRECTORIES)
    list(APPEND _arkode_example_directories ${_variant_dirs})
  endforeach()

  foreach(_example_dir ${_arkode_example_directories})
    file(RELATIVE_PATH _example_rel_dir ${CMAKE_CURRENT_SOURCE_DIR}
         ${_example_dir})
    sundials_install_example_cmake(
      arkode "${_example_dir}"
      "${SUNDIALS_EXAMPLES_INSTALL_PATH}/arkode/${_example_rel_dir}")
    sundials_add_test_install(arkode ${_example_rel_dir})
  endforeach()
endif()
