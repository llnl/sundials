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
foreach(example cvs_roberts_asai_klu cvs_roberts_klu)
  sundials_add_example(${example} SOURCES ${example}.c)
  target_link_libraries(
    ${example} PRIVATE ${CVODES_C_LIBS} sundials_nvecserial
                       sundials_sunlinsolklu ${EXE_EXTRA_LINK_LIBS})
endforeach()
sundials_add_example(
  cvs_roberts_fsa_klu SOURCES cvs_roberts_fsa_klu.c
  TEST_NAME cvs_roberts_fsa_klu_-sensi_stg1_t
  TEST_ARGS -sensi stg1 t
  LINK_LIBRARIES PRIVATE ${CVODES_C_LIBS} sundials_nvecserial
                 sundials_sunlinsolklu ${EXE_EXTRA_LINK_LIBS})
