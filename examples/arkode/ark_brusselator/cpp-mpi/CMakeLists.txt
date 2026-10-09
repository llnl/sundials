# -----------------------------------------------------------------
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
# -----------------------------------------------------------------
# ARKODE RAJA/MPI+X C++ Brusselator example

set(_source ark_brusselator1D_task_local_nls.cpp)
set(_common_files ark_brusselator1D.h plot_brusselator1D.py)
set(_common_libs sundials_arkode sundials_nvecmpiplusx RAJA
                 ${EXE_EXTRA_LINK_LIBS})

if(RAJA_BACKENDS MATCHES "HIP")
  list(APPEND _common_libs hip::device)
endif()
if((RAJA_BACKENDS MATCHES "TARGET_OPENMP") OR (RAJA_BACKENDS MATCHES "OPENMP"))
  list(APPEND _common_libs OpenMP::OpenMP_CXX)
endif()

# The serial NVECTOR variants have reference output for three command-line
# configurations.  The source is shared by all variants, so each backend gets
# a distinct target name.
sundials_add_example(
  ark_brusselator1D_task_local_nls.USE_SERIAL_NVEC ${_source}
  TEST_NAME ark_brusselator1D_task_local_nls.USE_SERIAL_NVEC_--monitor
  TEST_ARGS --monitor MPI_NPROCS 4
  ANSWER_FILE ark_brusselator1D_task_local_nls.USE_SERIAL_NVEC_--monitor.out
  EXAMPLE_TYPE develop INSTALL_FILES ${_common_files}
  ark_brusselator1D_task_local_nls.USE_SERIAL_NVEC_--monitor_--global-nls.out
  ark_brusselator1D_task_local_nls.USE_SERIAL_NVEC_--monitor_--explicit_--tf_1.out
  COMPILE_DEFINITIONS PRIVATE USE_SERIAL_NVEC
  LINK_LIBRARIES PRIVATE MPI::MPI_CXX ${_common_libs})

foreach(_args "--monitor --global-nls" "--monitor --explicit --tf 1")
  string(REPLACE " " "_" _suffix "${_args}")
  sundials_add_test(
    ark_brusselator1D_task_local_nls.USE_SERIAL_NVEC_${_suffix}
    ark_brusselator1D_task_local_nls.USE_SERIAL_NVEC TEST_ARGS ${_args}
    MPI_NPROCS 4 ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ark_brusselator1D_task_local_nls.USE_SERIAL_NVEC_${_suffix}.out
    EXAMPLE_TYPE develop)
endforeach()

# CUDA NVECTOR variants use the CUDA compiler for the shared C++ source.
if(SUNDIALS_ENABLE_NVECTOR_CUDA)
  set_source_files_properties(${_source} PROPERTIES LANGUAGE CUDA)
  set(_cuda_libs ${_common_libs} sundials_nveccuda)

  sundials_add_example(
    ark_brusselator1D_task_local_nls.USE_CUDA_NVEC ${_source}
    TEST_NAME ark_brusselator1D_task_local_nls.USE_CUDA_NVEC_--monitor
    TEST_ARGS --monitor MPI_NPROCS 4
    ANSWER_FILE ark_brusselator1D_task_local_nls.USE_CUDA_NVEC_--monitor.out
    EXAMPLE_TYPE develop INSTALL_FILES ${_common_files}
    ark_brusselator1D_task_local_nls.USE_CUDA_NVEC_--monitor_--global-nls.out
    COMPILE_DEFINITIONS PRIVATE USE_CUDA_NVEC
    LINK_LIBRARIES PRIVATE MPI::MPI_CXX ${_cuda_libs})

  sundials_add_test(
    ark_brusselator1D_task_local_nls.USE_CUDA_NVEC_--monitor_--global-nls
    ark_brusselator1D_task_local_nls.USE_CUDA_NVEC TEST_ARGS --monitor --global-nls
    MPI_NPROCS 4 ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ark_brusselator1D_task_local_nls.USE_CUDA_NVEC_--monitor_--global-nls.out
    EXAMPLE_TYPE develop)

  # These backends were built by the former CXX_parallel example when the
  # corresponding RAJA/NVECTOR support was enabled. They have no checked-in
  # answer file, so keep them as buildable, installed examples without tests.
  sundials_add_example(
    ark_brusselator1D_task_local_nls.USE_CUDAUVM_NVEC ${_source} NO_TEST
    INSTALL_FILES ${_common_files}
    COMPILE_DEFINITIONS PRIVATE USE_CUDAUVM_NVEC
    LINK_LIBRARIES PRIVATE MPI::MPI_CXX ${_cuda_libs})
endif()

if(SUNDIALS_ENABLE_NVECTOR_RAJA AND SUNDIALS_RAJA_BACKENDS MATCHES "CUDA")
  sundials_add_example(
    ark_brusselator1D_task_local_nls.USE_RAJA_NVEC ${_source} NO_TEST
    INSTALL_FILES ${_common_files}
    COMPILE_DEFINITIONS PRIVATE USE_RAJA_NVEC
    LINK_LIBRARIES PRIVATE MPI::MPI_CXX ${_common_libs} sundials_nvecraja)
endif()

if(SUNDIALS_ENABLE_NVECTOR_HIP AND RAJA_BACKENDS MATCHES "HIP")
  sundials_add_example(
    ark_brusselator1D_task_local_nls.USE_HIP_NVEC ${_source} NO_TEST
    INSTALL_FILES ${_common_files}
    COMPILE_DEFINITIONS PRIVATE USE_HIP_NVEC
    LINK_LIBRARIES PRIVATE MPI::MPI_CXX ${_common_libs} sundials_nvechip)
endif()

if(SUNDIALS_ENABLE_NVECTOR_OPENMPDEV AND RAJA_BACKENDS MATCHES "TARGET_OPENMP")
  sundials_add_example(
    ark_brusselator1D_task_local_nls.USE_OMPDEV_NVEC ${_source} NO_TEST
    INSTALL_FILES ${_common_files}
    COMPILE_DEFINITIONS PRIVATE USE_OMPDEV_NVEC
    LINK_LIBRARIES PRIVATE MPI::MPI_CXX ${_common_libs}
                   sundials_nvecopenmpdev)
endif()
