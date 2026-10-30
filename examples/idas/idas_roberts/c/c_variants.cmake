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
foreach(example idas_roberts_asai_dns idas_roberts_dns)
  set(example_options)
  if(example STREQUAL "idas_roberts_asai_dns")
    list(
      APPEND
      example_options
      EXAMPLE_TYPE
      exclude-single
      INSTALL_FILES
      idas_roberts_asai_dns_bkw1_stats.csv
      idas_roberts_asai_dns_fwd_stats.csv)
  else()
    list(APPEND example_options NO_TEST INSTALL_FILES idas_roberts_dns.out
         idas_roberts_dns_stats.csv)
  endif()
  sundials_add_example(${example} ${example}.c ${example_options})
  target_link_libraries(${example} PRIVATE ${IDAS_C_LIBS} sundials_nvecserial
                                           ${EXE_EXTRA_LINK_LIBS})
endforeach()
sundials_add_test(
  idas_roberts_dns idas_roberts_dns
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE idas_roberts_dns.out)

sundials_add_example(
  idas_roberts_fsa_dns idas_roberts_fsa_dns.c
  TEST_NAME idas_roberts_fsa_dns_-sensi_stg_t
  TEST_ARGS -sensi stg t
  EXAMPLE_TYPE exclude-single
  INSTALL_FILES idas_roberts_fsa_dns_stats_-sensi_stg_t.csv)
target_link_libraries(
  idas_roberts_fsa_dns PRIVATE ${IDAS_C_LIBS} sundials_nvecserial
                               ${EXE_EXTRA_LINK_LIBS})
