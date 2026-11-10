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
set(example_options TEST_ARGS 4)
sundials_add_example(
  idas_food_web_bnd_omp
  EXAMPLE_TYPE develop
  SOURCES idas_food_web_bnd_omp.c ${example_options}
  LINK_LIBRARIES ${IDAS_C_LIBS} sundials_nvecopenmp ${EXE_EXTRA_LINK_LIBS})
