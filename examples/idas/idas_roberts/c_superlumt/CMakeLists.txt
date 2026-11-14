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
set(example_options EXAMPLE_TYPE exclude-single)
set(answer_file idas_roberts_sps.out)
# Reference output was generated with 32-bit indices. Skip comparison with
# 64-bit indices until index-specific answers are added.
if(NOT SUNDIALS_INDEX_SIZE MATCHES "32")
  list(APPEND example_options NO_TEST INSTALL_FILES ${answer_file})
endif()
sundials_add_example(
  idas_roberts_sps
  SOURCES idas_roberts_sps.c ${example_options}
  LINK_LIBRARIES ${IDAS_C_LIBS} sundials_nvecserial sundials_sunlinsolsuperlumt
                 ${EXE_EXTRA_LINK_LIBS})
