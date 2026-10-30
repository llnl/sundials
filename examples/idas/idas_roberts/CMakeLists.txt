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
if(SUNDIALS_ENABLE_C_EXAMPLES)
  add_subdirectory(c_klu_roberts_asai_klu)
  add_subdirectory(c_klu_roberts_fsa_klu)
  add_subdirectory(c_klu_roberts_klu)
  add_subdirectory(c_roberts_asai_dns)
  add_subdirectory(c_roberts_dns)
  add_subdirectory(c_roberts_fsa_dns)
  add_subdirectory(c_superlumt)
  if(SUNDIALS_ENABLE_SUNLINSOL_KLU)
    add_subdirectory(c_klu)
  endif()
  if(SUNDIALS_ENABLE_SUNLINSOL_SUPERLUMT)
    add_subdirectory(c_superlumt)
  endif()
endif()
