# ---------------------------------------------------------------
# Programmer(s): David J. Gardner, Slaven Peles, and
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
# examples/cvode level CMakeLists.txt for SUNDIALS
# ---------------------------------------------------------------

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES README.md DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode)
endif()

# C examples
if(SUNDIALS_ENABLE_C_EXAMPLES)
  add_subdirectory(c)
endif()

# CXX examples
if(SUNDIALS_ENABLE_CXX_EXAMPLES)
  add_subdirectory(cpp)
  if(SUNDIALS_ENABLE_HIP)
    add_subdirectory(hip)
  endif()
  if(SUNDIALS_ENABLE_NVECTOR_SYCL)
    add_subdirectory(sycl)
  endif()
endif()

# Fortran examples
if(SUNDIALS_ENABLE_FORTRAN_EXAMPLES)
  add_subdirectory(fortran)
endif()

# CUDA examples
if(SUNDIALS_ENABLE_CUDA_EXAMPLES)

  if(SUNDIALS_ENABLE_CUDA AND CMAKE_CUDA_COMPILER)
    add_subdirectory(cuda)
  endif()

endif()
