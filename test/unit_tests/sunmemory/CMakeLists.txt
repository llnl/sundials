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

add_subdirectory(sys)

if(SUNDIALS_ENABLE_CUDA)
  add_subdirectory(cuda)
endif()

if(SUNDIALS_ENABLE_HIP)
  add_subdirectory(hip)
endif()

if(SUNDIALS_ENABLE_SYCL)
  add_subdirectory(sycl)
endif()
