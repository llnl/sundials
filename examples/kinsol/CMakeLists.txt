# ---------------------------------------------------------------
# Programmer(s): David J. Gardner and Slaven Peles @ LLNL
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
# examples/kinsol level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

set(KINSOL_C_LIBS sundials_kinsol)

# Examples are grouped by problem. Each problem directory contains variants
# grouped by implementation language and optional TPLs.
add_subdirectory(kin_analytic)
add_subdirectory(kin_bratu_2d)
add_subdirectory(kin_diagonal)
add_subdirectory(kin_em)
add_subdirectory(kin_ferraris_tronconi)
add_subdirectory(kin_food_web)
add_subdirectory(kin_heat_2d)
add_subdirectory(kin_krylov_demo)
add_subdirectory(kin_laplace)
add_subdirectory(kin_roberts)
add_subdirectory(kin_robot_kinematics)

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  get_property(
    _kinsol_problem_directories
    DIRECTORY
    PROPERTY SUBDIRECTORIES)
  set(_kinsol_example_directories)
  foreach(_problem_dir ${_kinsol_problem_directories})
    get_property(
      _variant_dirs
      DIRECTORY ${_problem_dir}
      PROPERTY SUBDIRECTORIES)
    list(APPEND _kinsol_example_directories ${_variant_dirs})
  endforeach()

  foreach(_example_dir ${_kinsol_example_directories})
    file(RELATIVE_PATH _example_rel_dir ${CMAKE_CURRENT_SOURCE_DIR}
         ${_example_dir})
    if(_example_rel_dir STREQUAL "kin_analytic/c")
      sundials_add_test_install(kinsol ${_example_rel_dir}
                                EXECUTABLE kin_analytic_fp)
    else()
      sundials_add_test_install(kinsol ${_example_rel_dir})
    endif()
  endforeach()
endif()
