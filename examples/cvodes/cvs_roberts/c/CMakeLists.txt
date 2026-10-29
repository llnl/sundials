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
  cvs_roberts_asai_dns cvs_roberts_asai_dns.c
  EXAMPLE_TYPE exclude-single
  INSTALL_FILES
    cvs_roberts_asai_dns_bkw1_stats.csv cvs_roberts_asai_dns_bkw2_stats.csv
    cvs_roberts_asai_dns_fwd_stats.csv
  LINK_LIBRARIES PRIVATE ${CVODES_C_LIBS} sundials_nvecserial
                 ${EXE_EXTRA_LINK_LIBS})
sundials_add_example(
  cvs_roberts_asai_dns_constraints cvs_roberts_asai_dns_constraints.c
  LINK_LIBRARIES PRIVATE ${CVODES_C_LIBS} sundials_nvecserial
                 ${EXE_EXTRA_LINK_LIBS})

sundials_add_example(
  cvs_roberts_fsa_dns cvs_roberts_fsa_dns.c
  TEST_NAME cvs_roberts_fsa_dns_-sensi_sim_t
  TEST_ARGS -sensi sim t
  EXAMPLE_TYPE exclude-single
  INSTALL_FILES
    cvs_roberts_fsa_dns_-sensi_stg1_t.out
    cvs_roberts_fsa_dns_stats_-sensi_sim_t.csv
    cvs_roberts_fsa_dns_stats_-sensi_stg1_t.csv
  LINK_LIBRARIES PRIVATE ${CVODES_C_LIBS} sundials_nvecserial
                 ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(
  cvs_roberts_fsa_dns_-sensi_stg1_t cvs_roberts_fsa_dns
  TEST_ARGS -sensi stg1 t
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE cvs_roberts_fsa_dns_-sensi_stg1_t.out
  EXAMPLE_TYPE exclude-single)

sundials_add_example(cvs_roberts_fsa_dns_switch cvs_roberts_fsa_dns_switch.c
                     EXAMPLE_TYPE exclude-single)
target_link_libraries(
  cvs_roberts_fsa_dns_switch PRIVATE ${CVODES_C_LIBS} sundials_nvecserial
                                     ${EXE_EXTRA_LINK_LIBS})
sundials_add_example(
  cvs_roberts_fsa_dns_constraints cvs_roberts_fsa_dns_constraints.c
  TEST_NAME cvs_roberts_fsa_dns_constraints_-sensi_stg1_t
  TEST_ARGS -sensi stg1 t
  EXAMPLE_TYPE exclude-single
  LINK_LIBRARIES PRIVATE ${CVODES_C_LIBS} sundials_nvecserial
                 ${EXE_EXTRA_LINK_LIBS})

sundials_add_example(
  cvs_roberts_dns cvs_roberts_dns.c
  INSTALL_FILES cvs_roberts_dns_stats.csv
  LINK_LIBRARIES PRIVATE ${CVODES_C_LIBS} sundials_nvecserial
                 ${EXE_EXTRA_LINK_LIBS})
foreach(example cvs_roberts_dns_constraints cvs_roberts_dns_uw)
  sundials_add_example(${example} ${example}.c)
  target_link_libraries(${example} PRIVATE ${CVODES_C_LIBS} sundials_nvecserial
                                           ${EXE_EXTRA_LINK_LIBS})
endforeach()
