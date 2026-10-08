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
foreach(_example ark_onewaycouple_mri ark_reaction_diffusion_mri
                 ark_twowaycouple_mri)
  set(_files)
  if(_example STREQUAL "ark_reaction_diffusion_mri")
    set(_files ark_reaction_diffusion_mri_fast_stats.csv
               ark_reaction_diffusion_mri_slow_stats.csv)
  endif()
  sundials_add_example(${_example} ${_example}.c EXAMPLE_TYPE develop
                        INSTALL_FILES ${_files}
                        LINK_LIBRARIES PRIVATE sundials_arkode
                        sundials_nvecserial ${EXE_EXTRA_LINK_LIBS})
endforeach()
