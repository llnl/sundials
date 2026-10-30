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
sundials_add_example(
  ida_food_web_bnd_omp SOURCES ida_food_web_bnd_omp.c
  TEST_ARGS 4
  EXAMPLE_TYPE develop
  LINK_LIBRARIES PRIVATE ${IDA_C_LIBS} sundials_nvecopenmp
                 ${EXE_EXTRA_LINK_LIBS})
sundials_add_example(
  ida_food_web_kry_omp SOURCES ida_food_web_kry_omp.c NO_TEST
  INSTALL_FILES ida_food_web_kry_omp.out
  LINK_LIBRARIES PRIVATE ${IDA_C_LIBS} sundials_nvecopenmp
                 ${EXE_EXTRA_LINK_LIBS})
