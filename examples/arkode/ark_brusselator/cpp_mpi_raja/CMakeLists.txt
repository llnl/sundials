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
# Tests are registered separately below for their argument-specific answers.
set(_brusselator_task_local_source
    "${CMAKE_CURRENT_SOURCE_DIR}/ark_brusselator_1d_task_local_nls.cpp")
set(_brusselator_task_local_header
    "${CMAKE_CURRENT_SOURCE_DIR}/ark_brusselator_1d.h")
set(_common_libs sundials_arkode sundials_nvecmpiplusx RAJA
                 ${EXE_EXTRA_LINK_LIBS})
if(RAJA_BACKENDS MATCHES "HIP")
  list(APPEND _common_libs hip::device)
endif()
if((RAJA_BACKENDS MATCHES "TARGET_OPENMP") OR (RAJA_BACKENDS MATCHES "OPENMP"))
  list(APPEND _common_libs OpenMP::OpenMP_CXX)
endif()
sundials_add_example(
  ark_brusselator_1d_task_local_nls.USE_SERIAL_NVEC
  SOURCES ${_brusselator_task_local_source}
  TEST_NAME ark_brusselator_1d_task_local_nls.USE_SERIAL_NVEC_--monitor
  TEST_ARGS --monitor
  MPI_NPROCS 4
  ANSWER_FILE ark_brusselator_1d_task_local_nls.USE_SERIAL_NVEC_--monitor.out
  EXAMPLE_TYPE develop
  INSTALL_FILES
    ${_brusselator_task_local_header}
    plot_brusselator_1d.py
    ark_brusselator_1d_task_local_nls.USE_SERIAL_NVEC_--monitor_--global-nls.out
    ark_brusselator_1d_task_local_nls.USE_SERIAL_NVEC_--monitor_--explicit_--tf_1.out
  COMPILE_DEFINITIONS PRIVATE USE_SERIAL_NVEC
  LINK_LIBRARIES PRIVATE MPI::MPI_CXX ${_common_libs})
foreach(_args "--monitor --global-nls" "--monitor --explicit --tf 1")
  string(REPLACE " " "_" _suffix "${_args}")
  sundials_add_test(
    ark_brusselator_1d_task_local_nls.USE_SERIAL_NVEC_${_suffix}
    ark_brusselator_1d_task_local_nls.USE_SERIAL_NVEC
    TEST_ARGS ${_args}
    MPI_NPROCS 4
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ark_brusselator_1d_task_local_nls.USE_SERIAL_NVEC_${_suffix}.out
    EXAMPLE_TYPE develop)
endforeach()
if(SUNDIALS_ENABLE_NVECTOR_CUDA)
  set_source_files_properties(${_brusselator_task_local_source}
                              PROPERTIES LANGUAGE CUDA)
  set(_cuda_libs ${_common_libs} sundials_nveccuda)
endif()

set(_common_libs sundials_arkode sundials_nvecmpiplusx RAJA
                 ${EXE_EXTRA_LINK_LIBS})
if(RAJA_BACKENDS MATCHES "HIP")
  list(APPEND _common_libs hip::device)
endif()
if((RAJA_BACKENDS MATCHES "TARGET_OPENMP") OR (RAJA_BACKENDS MATCHES "OPENMP"))
  list(APPEND _common_libs OpenMP::OpenMP_CXX)
endif()
if(SUNDIALS_ENABLE_NVECTOR_CUDA)
  set_source_files_properties(${_brusselator_task_local_source}
                              PROPERTIES LANGUAGE CUDA)
  set(_cuda_libs ${_common_libs} sundials_nveccuda)
  sundials_add_example(
    ark_brusselator_1d_task_local_nls.USE_CUDA_NVEC
    SOURCES ${_brusselator_task_local_source}
    TEST_NAME ark_brusselator_1d_task_local_nls.USE_CUDA_NVEC_--monitor
    TEST_ARGS --monitor
    MPI_NPROCS 4
    ANSWER_FILE ark_brusselator_1d_task_local_nls.USE_CUDA_NVEC_--monitor.out
    EXAMPLE_TYPE develop
    INSTALL_FILES
      ${_brusselator_task_local_header} plot_brusselator_1d.py
      ark_brusselator_1d_task_local_nls.USE_CUDA_NVEC_--monitor_--global-nls.out
    COMPILE_DEFINITIONS PRIVATE USE_CUDA_NVEC
    LINK_LIBRARIES PRIVATE MPI::MPI_CXX ${_cuda_libs})
  sundials_add_test(
    ark_brusselator_1d_task_local_nls.USE_CUDA_NVEC_--monitor_--global-nls
    ark_brusselator_1d_task_local_nls.USE_CUDA_NVEC
    TEST_ARGS --monitor --global-nls
    MPI_NPROCS 4
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE
      ark_brusselator_1d_task_local_nls.USE_CUDA_NVEC_--monitor_--global-nls.out
    EXAMPLE_TYPE develop)
endif()

set(_common_libs sundials_arkode sundials_nvecmpiplusx RAJA
                 ${EXE_EXTRA_LINK_LIBS})
if(RAJA_BACKENDS MATCHES "HIP")
  list(APPEND _common_libs hip::device)
endif()
if((RAJA_BACKENDS MATCHES "TARGET_OPENMP") OR (RAJA_BACKENDS MATCHES "OPENMP"))
  list(APPEND _common_libs OpenMP::OpenMP_CXX)
endif()
if(SUNDIALS_ENABLE_NVECTOR_CUDA)
  set_source_files_properties(${_brusselator_task_local_source}
                              PROPERTIES LANGUAGE CUDA)
  set(_cuda_libs ${_common_libs} sundials_nveccuda)
  sundials_add_example(
    ark_brusselator_1d_task_local_nls.USE_CUDAUVM_NVEC
    SOURCES ${_brusselator_task_local_source}
    NO_TEST
    INSTALL_FILES ${_brusselator_task_local_header} plot_brusselator_1d.py
    COMPILE_DEFINITIONS PRIVATE USE_CUDAUVM_NVEC
    LINK_LIBRARIES PRIVATE MPI::MPI_CXX ${_cuda_libs})
endif()

set(_common_libs sundials_arkode sundials_nvecmpiplusx RAJA
                 ${EXE_EXTRA_LINK_LIBS})
if(RAJA_BACKENDS MATCHES "HIP")
  list(APPEND _common_libs hip::device)
endif()
if((RAJA_BACKENDS MATCHES "TARGET_OPENMP") OR (RAJA_BACKENDS MATCHES "OPENMP"))
  list(APPEND _common_libs OpenMP::OpenMP_CXX)
endif()
if(SUNDIALS_ENABLE_NVECTOR_CUDA)
  set_source_files_properties(${_brusselator_task_local_source}
                              PROPERTIES LANGUAGE CUDA)
  set(_cuda_libs ${_common_libs} sundials_nveccuda)
endif()
if(SUNDIALS_ENABLE_NVECTOR_HIP AND RAJA_BACKENDS MATCHES "HIP")
  sundials_add_example(
    ark_brusselator_1d_task_local_nls.USE_HIP_NVEC
    SOURCES ${_brusselator_task_local_source}
    NO_TEST
    INSTALL_FILES ${_brusselator_task_local_header} plot_brusselator_1d.py
    COMPILE_DEFINITIONS PRIVATE USE_HIP_NVEC
    LINK_LIBRARIES PRIVATE MPI::MPI_CXX ${_common_libs} sundials_nvechip)
endif()

set(_common_libs sundials_arkode sundials_nvecmpiplusx RAJA
                 ${EXE_EXTRA_LINK_LIBS})
if(RAJA_BACKENDS MATCHES "HIP")
  list(APPEND _common_libs hip::device)
endif()
if((RAJA_BACKENDS MATCHES "TARGET_OPENMP") OR (RAJA_BACKENDS MATCHES "OPENMP"))
  list(APPEND _common_libs OpenMP::OpenMP_CXX)
endif()
if(SUNDIALS_ENABLE_NVECTOR_CUDA)
  set_source_files_properties(${_brusselator_task_local_source}
                              PROPERTIES LANGUAGE CUDA)
  set(_cuda_libs ${_common_libs} sundials_nveccuda)
endif()
if(SUNDIALS_ENABLE_NVECTOR_OPENMPDEV AND RAJA_BACKENDS MATCHES "TARGET_OPENMP")
  sundials_add_example(
    ark_brusselator_1d_task_local_nls.USE_OMPDEV_NVEC
    SOURCES ${_brusselator_task_local_source}
    NO_TEST
    INSTALL_FILES ${_brusselator_task_local_header} plot_brusselator_1d.py
    COMPILE_DEFINITIONS PRIVATE USE_OMPDEV_NVEC
    LINK_LIBRARIES PRIVATE MPI::MPI_CXX ${_common_libs} sundials_nvecopenmpdev)
endif()

set(_common_libs sundials_arkode sundials_nvecmpiplusx RAJA
                 ${EXE_EXTRA_LINK_LIBS})
if(RAJA_BACKENDS MATCHES "HIP")
  list(APPEND _common_libs hip::device)
endif()
if((RAJA_BACKENDS MATCHES "TARGET_OPENMP") OR (RAJA_BACKENDS MATCHES "OPENMP"))
  list(APPEND _common_libs OpenMP::OpenMP_CXX)
endif()
if(SUNDIALS_ENABLE_NVECTOR_CUDA)
  set_source_files_properties(${_brusselator_task_local_source}
                              PROPERTIES LANGUAGE CUDA)
  set(_cuda_libs ${_common_libs} sundials_nveccuda)
endif()
if(SUNDIALS_ENABLE_NVECTOR_RAJA AND SUNDIALS_RAJA_BACKENDS MATCHES "CUDA")
  sundials_add_example(
    ark_brusselator_1d_task_local_nls.USE_RAJA_NVEC
    SOURCES ${_brusselator_task_local_source}
    NO_TEST
    INSTALL_FILES ${_brusselator_task_local_header} plot_brusselator_1d.py
    COMPILE_DEFINITIONS PRIVATE USE_RAJA_NVEC
    LINK_LIBRARIES PRIVATE MPI::MPI_CXX ${_common_libs} sundials_nvecraja)
endif()
