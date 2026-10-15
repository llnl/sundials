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
foreach(example idasRoberts_ASAi_klu idasRoberts_klu)
  sundials_add_example(${example} ${example}.c)
  target_link_libraries(
    ${example} PRIVATE ${IDAS_C_LIBS} sundials_nvecserial
                       sundials_sunlinsolklu ${EXE_EXTRA_LINK_LIBS})
endforeach()

sundials_add_example(idasRoberts_FSA_klu idasRoberts_FSA_klu.c
                     TEST_NAME idasRoberts_FSA_klu_-sensi_stg_t
                     TEST_ARGS -sensi stg t)
target_link_libraries(
  idasRoberts_FSA_klu PRIVATE ${IDAS_C_LIBS} sundials_nvecserial
                              sundials_sunlinsolklu ${EXE_EXTRA_LINK_LIBS})
