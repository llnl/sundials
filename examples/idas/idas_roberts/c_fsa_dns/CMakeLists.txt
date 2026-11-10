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
# Preserve develop's single-precision test exclusion for these cases.
sundials_add_example(
  idas_roberts_fsa_dns
  SOURCES idas_roberts_fsa_dns.c
  TEST_NAME idas_roberts_fsa_dns_-sensi_stg_t
  TEST_ARGS -sensi stg t
  EXAMPLE_TYPE exclude-single
  INSTALL_FILES idas_roberts_fsa_dns_stats_-sensi_stg_t.csv
  LINK_LIBRARIES PRIVATE ${IDAS_C_LIBS} sundials_nvecserial
                 ${EXE_EXTRA_LINK_LIBS})
