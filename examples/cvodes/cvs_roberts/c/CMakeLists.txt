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
  cvsRoberts_ASAi_dns cvsRoberts_ASAi_dns.c
  EXAMPLE_TYPE exclude-single
  INSTALL_FILES cvsRoberts_ASAi_dns_bkw1_stats.csv
                cvsRoberts_ASAi_dns_bkw2_stats.csv
                cvsRoberts_ASAi_dns_fwd_stats.csv
  LINK_LIBRARIES PRIVATE ${CVODES_C_LIBS} sundials_nvecserial
                 ${EXE_EXTRA_LINK_LIBS})
sundials_add_example(
  cvsRoberts_ASAi_dns_constraints cvsRoberts_ASAi_dns_constraints.c
  LINK_LIBRARIES PRIVATE ${CVODES_C_LIBS} sundials_nvecserial
                 ${EXE_EXTRA_LINK_LIBS})

sundials_add_example(
  cvsRoberts_FSA_dns cvsRoberts_FSA_dns.c
  TEST_NAME cvsRoberts_FSA_dns_-sensi_sim_t
  TEST_ARGS -sensi sim t
  EXAMPLE_TYPE exclude-single
  INSTALL_FILES cvsRoberts_FSA_dns_-sensi_stg1_t.out
                cvsRoberts_FSA_dns_stats_-sensi_sim_t.csv
                cvsRoberts_FSA_dns_stats_-sensi_stg1_t.csv
  LINK_LIBRARIES PRIVATE ${CVODES_C_LIBS} sundials_nvecserial
                 ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(
  cvsRoberts_FSA_dns_-sensi_stg1_t cvsRoberts_FSA_dns
  TEST_ARGS -sensi stg1 t
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE cvsRoberts_FSA_dns_-sensi_stg1_t.out
  EXAMPLE_TYPE exclude-single)

sundials_add_example(cvsRoberts_FSA_dns_Switch cvsRoberts_FSA_dns_Switch.c
                     EXAMPLE_TYPE exclude-single)
target_link_libraries(cvsRoberts_FSA_dns_Switch PRIVATE ${CVODES_C_LIBS}
                                                       sundials_nvecserial
                                                       ${EXE_EXTRA_LINK_LIBS})
sundials_add_example(
  cvsRoberts_FSA_dns_constraints cvsRoberts_FSA_dns_constraints.c
  TEST_NAME cvsRoberts_FSA_dns_constraints_-sensi_stg1_t
  TEST_ARGS -sensi stg1 t
  EXAMPLE_TYPE exclude-single
  LINK_LIBRARIES PRIVATE ${CVODES_C_LIBS} sundials_nvecserial
                 ${EXE_EXTRA_LINK_LIBS})

sundials_add_example(
  cvsRoberts_dns cvsRoberts_dns.c
  INSTALL_FILES cvsRoberts_dns_stats.csv
  LINK_LIBRARIES PRIVATE ${CVODES_C_LIBS} sundials_nvecserial
                 ${EXE_EXTRA_LINK_LIBS})
foreach(example cvsRoberts_dns_constraints cvsRoberts_dns_uw)
  sundials_add_example(${example} ${example}.c)
  target_link_libraries(${example} PRIVATE ${CVODES_C_LIBS}
                                           sundials_nvecserial
                                           ${EXE_EXTRA_LINK_LIBS})
endforeach()
