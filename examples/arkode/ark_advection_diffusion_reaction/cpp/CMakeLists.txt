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
sundials_add_example(ark_adr1d_extsts ark_adr1d_extsts.cpp
                      EXAMPLE_TYPE develop INSTALL_FILES ark_adr1d_extsts.hpp
                      INCLUDE_DIRECTORIES PRIVATE ${PROJECT_SOURCE_DIR}/examples/utilities
                      LINK_LIBRARIES PRIVATE sundials_arkode sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})

if(SUNDIALS_ENABLE_CVODE)
  sundials_add_example(ark_advection_diffusion_reaction
                        ark_advection_diffusion_reaction.cpp
                        EXAMPLE_TYPE develop
                        INSTALL_FILES ark_advection_diffusion_reaction.hpp
                        INCLUDE_DIRECTORIES PRIVATE ${PROJECT_SOURCE_DIR}/examples/utilities
                        LINK_LIBRARIES PRIVATE sundials_arkode sundials_cvode
                        sundials_nvecserial ${EXE_EXTRA_LINK_LIBS})
endif()
