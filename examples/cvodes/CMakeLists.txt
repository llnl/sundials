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
# examples/cvodes level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

set(CVODES_C_LIBS sundials_cvodes)

# Examples are grouped by problem. Each problem directory contains variants
# grouped by implementation language and optional TPLs.
add_subdirectory(cvs_adv_diff)
add_subdirectory(cvs_analytic)
add_subdirectory(cvs_atmospheric_dispersion)
add_subdirectory(cvs_direct_demo)
add_subdirectory(cvs_diurnal)
add_subdirectory(cvs_food_web)
add_subdirectory(cvs_hessian)
add_subdirectory(cvs_krylov_demo)
add_subdirectory(cvs_lotka_volterra)
add_subdirectory(cvs_particle)
add_subdirectory(cvs_pendulum)
add_subdirectory(cvs_roberts)

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  get_property(
    _cvodes_problem_directories
    DIRECTORY
    PROPERTY SUBDIRECTORIES)
  set(_cvodes_example_directories)
  foreach(_problem_dir ${_cvodes_problem_directories})
    get_property(
      _variant_dirs
      DIRECTORY ${_problem_dir}
      PROPERTY SUBDIRECTORIES)
    list(APPEND _cvodes_example_directories ${_variant_dirs})
  endforeach()

  foreach(_example_dir ${_cvodes_example_directories})
    file(RELATIVE_PATH _example_rel_dir ${CMAKE_CURRENT_SOURCE_DIR}
         ${_example_dir})
    if(_example_rel_dir STREQUAL "cvs_roberts/c")
      sundials_add_test_install(cvodes ${_example_rel_dir}
                                EXECUTABLE cvs_roberts_dns)
    else()
      sundials_add_test_install(cvodes ${_example_rel_dir})
    endif()
  endforeach()
endif()
