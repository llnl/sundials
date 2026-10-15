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
foreach(example idasRoberts_ASAi_dns idasRoberts_dns)
  set(example_options)
  if(example STREQUAL "idasRoberts_ASAi_dns")
    list(APPEND example_options EXAMPLE_TYPE exclude-single INSTALL_FILES
         idasRoberts_ASAi_dns_bkw1_stats.csv
         idasRoberts_ASAi_dns_fwd_stats.csv)
  else()
    list(APPEND example_options INSTALL_FILES idasRoberts_dns_stats.csv)
  endif()
  sundials_add_example(${example} ${example}.c ${example_options})
  target_link_libraries(${example} PRIVATE ${IDAS_C_LIBS} sundials_nvecserial
                                           ${EXE_EXTRA_LINK_LIBS})
endforeach()

sundials_add_example(
  idasRoberts_FSA_dns idasRoberts_FSA_dns.c
  TEST_NAME idasRoberts_FSA_dns_-sensi_stg_t
  TEST_ARGS -sensi stg t
  EXAMPLE_TYPE exclude-single
  INSTALL_FILES idasRoberts_FSA_dns_stats_-sensi_stg_t.csv)
target_link_libraries(idasRoberts_FSA_dns PRIVATE ${IDAS_C_LIBS}
                                                 sundials_nvecserial
                                                 ${EXE_EXTRA_LINK_LIBS})
