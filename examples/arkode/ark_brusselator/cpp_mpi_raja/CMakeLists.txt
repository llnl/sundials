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
set(_common_libs sundials_arkode sundials_nvecmpiplusx RAJA
                 ${EXE_EXTRA_LINK_LIBS})
if(RAJA_BACKENDS MATCHES "HIP")
  list(APPEND _common_libs hip::device)
endif()
if((RAJA_BACKENDS MATCHES "TARGET_OPENMP") OR (RAJA_BACKENDS MATCHES "OPENMP"))
  list(APPEND _common_libs OpenMP::OpenMP_CXX)
endif()
if(SUNDIALS_ENABLE_NVECTOR_CUDA)
  set_source_files_properties(ark_brusselator_1d_task_local_nls.cpp
                              PROPERTIES LANGUAGE CUDA)
  set(_cuda_libs ${_common_libs} sundials_nveccuda)
endif()
if(SUNDIALS_ENABLE_NVECTOR_RAJA AND SUNDIALS_RAJA_BACKENDS MATCHES "CUDA")
  sundials_add_example(
    ark_brusselator_1d_task_local_nls.USE_RAJA_NVEC
    SOURCES ark_brusselator_1d_task_local_nls.cpp
    NO_TEST
    INSTALL_FILES ark_brusselator_1d.h plot_brusselator_1d.py
    COMPILE_DEFINITIONS PRIVATE USE_RAJA_NVEC
    LINK_LIBRARIES PRIVATE MPI::MPI_CXX ${_common_libs} sundials_nvecraja)
endif()
