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
sundials_add_example(cvsParticle_dns cvsParticle_dns.c
                     INSTALL_FILES plot_cvsParticle.py)
target_link_libraries(cvsParticle_dns PRIVATE ${CVODES_C_LIBS}
                                              sundials_nvecserial
                                              ${EXE_EXTRA_LINK_LIBS})
