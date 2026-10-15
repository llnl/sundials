# ---------------------------------------------------------------
# Programmer(s): David J. Gardner, Slaven Peles, and
#                Cody Balos @ LLNL
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
# examples/ida level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

set(IDA_C_LIBS sundials_ida)

add_subdirectory(ida_analytic)
add_subdirectory(ida_food_web)
add_subdirectory(ida_heat_2d)
add_subdirectory(ida_krylov_demo)
add_subdirectory(ida_roberts)
add_subdirectory(ida_slider_crank)

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  include(SundialsInstallExampleCMake)

  get_property(
    _ida_problem_directories
    DIRECTORY
    PROPERTY SUBDIRECTORIES)
  set(_ida_example_directories)
  foreach(_problem_dir ${_ida_problem_directories})
    get_property(
      _variant_dirs
      DIRECTORY ${_problem_dir}
      PROPERTY SUBDIRECTORIES)
    list(APPEND _ida_example_directories ${_variant_dirs})
  endforeach()

  foreach(_example_dir ${_ida_example_directories})
    file(RELATIVE_PATH _example_rel_dir ${CMAKE_CURRENT_SOURCE_DIR}
         ${_example_dir})
    sundials_install_example_cmake(
      ida "${_example_dir}"
      "${SUNDIALS_EXAMPLES_INSTALL_PATH}/ida/${_example_rel_dir}")
    sundials_add_test_install(ida ${_example_rel_dir})
  endforeach()
endif()
