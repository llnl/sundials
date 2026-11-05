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
set(example_options)
list(APPEND example_options NO_TEST INSTALL_FILES idas_roberts_dns.out
     idas_roberts_dns_stats.csv)
sundials_add_example(
  idas_roberts_dns
  SOURCES idas_roberts_dns.c ${example_options}
  LINK_LIBRARIES PRIVATE ${IDAS_C_LIBS} sundials_nvecserial
                 ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(
  idas_roberts_dns idas_roberts_dns
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE idas_roberts_dns.out)
