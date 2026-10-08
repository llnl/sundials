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
sundials_add_example(ark_bruss1D_FEM_klu_f2003 ark_bruss1D_FEM_klu_f2003.f90
                      EXAMPLE_TYPE develop
                      LINK_LIBRARIES PRIVATE sundials_farkode_mod sundials_arkode
                      sundials_sunlinsolklu sundials_fsunlinsolklu_mod
                      ${EXE_EXTRA_LINK_LIBS})
