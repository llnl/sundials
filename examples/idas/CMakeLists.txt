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
# examples/idas level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

set(IDAS_C_LIBS sundials_idas)

# Examples are grouped by problem. Each problem directory contains variants
# grouped by implementation language and optional TPLs.
add_subdirectory(idas_akzo_nobel)
add_subdirectory(idas_analytic)
add_subdirectory(idas_brusselator)
add_subdirectory(idas_food_web)
add_subdirectory(idas_heat_2d)
add_subdirectory(idas_hessian)
add_subdirectory(idas_krylov_demo)
add_subdirectory(idas_roberts)
add_subdirectory(idas_slider_crank)

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  include(SundialsInstallExampleCMake)

  get_property(
    _idas_problem_directories
    DIRECTORY
    PROPERTY SUBDIRECTORIES)
  set(_idas_example_directories)
  foreach(_problem_dir ${_idas_problem_directories})
    get_property(
      _variant_dirs
      DIRECTORY ${_problem_dir}
      PROPERTY SUBDIRECTORIES)
    list(APPEND _idas_example_directories ${_variant_dirs})
  endforeach()

  foreach(_example_dir ${_idas_example_directories})
    file(RELATIVE_PATH _example_rel_dir ${CMAKE_CURRENT_SOURCE_DIR}
         ${_example_dir})
    sundials_install_example_cmake(
      idas "${_example_dir}"
      "${SUNDIALS_EXAMPLES_INSTALL_PATH}/idas/${_example_rel_dir}")

    if(_example_rel_dir STREQUAL "idas_roberts/c")
      sundials_add_test_install(idas ${_example_rel_dir}
                                EXECUTABLE idasRoberts_dns)
    else()
      sundials_add_test_install(idas ${_example_rel_dir})
    endif()
  endforeach()
endif()
