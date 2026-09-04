# ---------------------------------------------------------------
# Programmer(s): Slaven Peles and Cody J. Balos @ LLNL
# ---------------------------------------------------------------
# Based on CMakeLists.txt in cvode/serial by Radu Serban @ LLNL
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
# CMakeLists.txt file for CVODE CUDA examples
# ---------------------------------------------------------------

# Example lists are tuples "name\;args\;type" where the type is 'develop' for
# examples excluded from 'make test' in releases

# Examples using SUNDIALS linear solvers
set(CVODE_examples
    "cvAdvDiff_kry_cuda\;\;develop" "cvAdvDiff_kry_cuda_managed\;\;develop"
    "cvAdvDiff_diag_cuda\;0 0\;develop" "cvAdvDiff_diag_cuda\;0 1\;develop"
    "cvAdvDiff_diag_cuda\;1 1\;develop")

# Examples using cuSolverSP linear solvers
set(CVODE_examples_cusolver "cvRoberts_block_cusolversp_batchqr\;\;develop")

if(SUNDIALS_INDEX_SIZE MATCHES "32")
  set(all_examples "${CVODE_examples_cusolver}")
  # Specify libraries to link against
  set(SUNDIALS_LIBS sundials_sunlinsolcusolversp)
endif()

list(APPEND all_examples "${CVODE_examples}")

# Specify libraries to link against
set(SUNDIALS_LIBS ${SUNDIALS_LIBS} sundials_cvode sundials_nveccuda)

if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
  set(SUNDIALS_LIBS ${SUNDIALS_LIBS} sundials_cvode_fused_cuda)
endif()

# Add source directory to include directories
include_directories(.)

# Build each example in its own directory. The example lists remain here for
# the installed CMake/Makefile templates below.
add_subdirectory(cvAdvDiff_kry_cuda)
add_subdirectory(cvAdvDiff_kry_cuda_managed)
add_subdirectory(cvAdvDiff_diag_cuda)
if(SUNDIALS_INDEX_SIZE MATCHES "32")
  add_subdirectory(cvRoberts_block_cusolversp_batchqr)
endif()

# create Makefile and CMakeLists.txt for examples
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)

  # Install the README file
  install(FILES README DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cuda)

  # Prepare substitution variables for Makefile and/or CMakeLists templates
  set(SOLVER "CVODE")
  set(SOLVER_LIB "sundials_cvode")
  set(NVECTOR_LIB "sundials_nveccuda")
  if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
    set(LIBS "-lsundials_cvode_fused_cuda ${LIBS}")
  endif()

  set(CVODE_examples_install
      "cvAdvDiff_kry_cuda/cvAdvDiff_kry_cuda\;\;develop"
      "cvAdvDiff_kry_cuda_managed/cvAdvDiff_kry_cuda_managed\;\;develop"
      "cvAdvDiff_diag_cuda/cvAdvDiff_diag_cuda\;0 0\;develop"
      "cvAdvDiff_diag_cuda/cvAdvDiff_diag_cuda\;0 1\;develop"
      "cvAdvDiff_diag_cuda/cvAdvDiff_diag_cuda\;1 1\;develop")
  examples2string(CVODE_examples_install EXAMPLES)

  if(SUNDIALS_INDEX_SIZE MATCHES "32")
    set(SUNLS_LIB "sundials_sunlinsolcusolversp")
    set(SUNMAT_LIB "sundials_sunmatrixcusparse")
    set(CVODE_examples_cusolver_install
        "cvRoberts_block_cusolversp_batchqr/cvRoberts_block_cusolversp_batchqr\;\;develop")
    examples2string(CVODE_examples_cusolver_install EXAMPLES_CUSOLVER)
  endif()

  # Regardless of the platform we're on, we will generate and install
  # CMakeLists.txt file for building the examples. This file  can then be used
  # as a template for the user's own programs.

  # generate CMakelists.txt in the binary directory
  configure_file(
    ${PROJECT_SOURCE_DIR}/examples/templates/cmakelists_serial_CUDA_ex.in
    ${PROJECT_BINARY_DIR}/examples/cvode/cuda/CMakeLists.txt @ONLY)

  # install CMakelists.txt
  install(FILES ${PROJECT_BINARY_DIR}/examples/cvode/cuda/CMakeLists.txt
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cuda)

  # On UNIX-type platforms, we also  generate and install a makefile for
  # building the examples. This makefile can then be used as a template for the
  # user's own programs.

  if(UNIX)
    # generate Makefile and place it in the binary dir
    configure_file(
      ${PROJECT_SOURCE_DIR}/examples/templates/makefile_serial_CUDA_ex.in
      ${PROJECT_BINARY_DIR}/examples/cvode/cuda/Makefile_ex @ONLY)
    # install the configured Makefile_ex as Makefile
    install(
      FILES ${PROJECT_BINARY_DIR}/examples/cvode/cuda/Makefile_ex
      DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cuda
      RENAME Makefile)
  endif()

  # add test_install target
  sundials_add_test_install(cvode cuda)

endif()
