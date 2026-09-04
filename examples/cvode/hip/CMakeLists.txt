# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# CMakeLists.txt file for CVODE HIP examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is 'develop' for
# examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(all_examples
    "cvAdvDiff_kry_hip.cpp\;\;develop" "cvAdvDiff_diag_hip.cpp\;0 0\;develop"
    "cvAdvDiff_diag_hip.cpp\;0 1\;develop"
    "cvAdvDiff_diag_hip.cpp\;1 1\;develop")

# Specify libraries to link against
set(SUNDIALS_LIBS sundials_cvode sundials_nvechip)

if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
  set(SUNDIALS_LIBS ${SUNDIALS_LIBS} sundials_cvode_fused_hip)
endif()

# Add source directory to include directories
include_directories(.)

# Build each example in its own directory. The example list remains here for
# the installed CMake template below.
add_subdirectory(cvAdvDiff_kry_hip)
add_subdirectory(cvAdvDiff_diag_hip)

# Install the targets
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)

  install(FILES README DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/hip)

  if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
    set(_fused_kernels_if_enabled cvode_fused_hip)
  else()
    set(_fused_kernels_if_enabled)
  endif()

  # Sources and answer files are installed by the leaf directories. Generate
  # the category-level installed CMake template here without asking the
  # install helper to look for the now-nested source files.
  set(SOLVER CVODE)
  set(EXAMPLES_DEPENDENCIES "")
  set(EXTRA_INCLUDES "")
  set(LIBS "")
  set(all_examples_install
      "cvAdvDiff_kry_hip/cvAdvDiff_kry_hip.cpp\;\;develop"
      "cvAdvDiff_diag_hip/cvAdvDiff_diag_hip.cpp\;0 0\;develop"
      "cvAdvDiff_diag_hip/cvAdvDiff_diag_hip.cpp\;0 1\;develop"
      "cvAdvDiff_diag_hip/cvAdvDiff_diag_hip.cpp\;1 1\;develop")
  examples2string(all_examples_install EXAMPLES)
  set(_components cvode nvechip ${_fused_kernels_if_enabled})
  set(EXAMPLES_CMAKE_COMPONENTS "")
  list2string(_components EXAMPLES_CMAKE_COMPONENTS)
  set(_targets "")
  foreach(_target ${_components})
    list(APPEND _targets SUNDIALS::${_target})
  endforeach()
  set(EXAMPLES_CMAKE_TARGETS "")
  list2string(_targets EXAMPLES_CMAKE_TARGETS)
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_HIP_ex.in
    ${PROJECT_BINARY_DIR}/examples/cvode/hip/CMakeLists.txt @ONLY)
  install(FILES ${PROJECT_BINARY_DIR}/examples/cvode/hip/CMakeLists.txt
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/hip)
  sundials_add_test_install(cvode hip)

endif()
