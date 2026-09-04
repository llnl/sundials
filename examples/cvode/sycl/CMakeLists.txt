# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
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
# CMakeLists.txt file for CVODE SYCL examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is 'develop' for
# examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(CVODE_examples "cvAdvDiff_kry_sycl.cpp\;\;develop")

# Build each example in its own directory. The example list remains here for
# the installed CMake template below.
add_subdirectory(cvAdvDiff_kry_sycl)

# Install the targets
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)

  install(FILES README DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/sycl)

  # Sources and answer files are installed by the leaf directory. Generate the
  # category-level installed CMake template without looking for flat sources.
  set(SOLVER CVODE)
  set(EXAMPLES_DEPENDENCIES "")
  set(EXTRA_INCLUDES "")
  set(EXAMPLES_FIND_PACKAGE "")
  set(EXAMPLES_EXTRA_LIBS "")
  set(CVODE_examples_install
      "cvAdvDiff_kry_sycl/cvAdvDiff_kry_sycl.cpp\;\;develop")
  examples2string(CVODE_examples_install EXAMPLES)
  set(_components cvode nvecsycl)
  set(EXAMPLES_CMAKE_COMPONENTS "")
  list2string(_components EXAMPLES_CMAKE_COMPONENTS)
  set(_targets SUNDIALS::cvode SUNDIALS::nvecsycl)
  set(EXAMPLES_CMAKE_TARGETS "")
  list2string(_targets EXAMPLES_CMAKE_TARGETS)
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_CXX_ex.in
    ${PROJECT_BINARY_DIR}/examples/cvode/sycl/CMakeLists.txt @ONLY)
  install(FILES ${PROJECT_BINARY_DIR}/examples/cvode/sycl/CMakeLists.txt
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/sycl)
  sundials_add_test_install(cvode sycl)

endif()
