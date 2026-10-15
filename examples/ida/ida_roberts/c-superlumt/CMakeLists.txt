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
if(NOT SUNDIALS_INDEX_SIZE MATCHES "32")
  list(APPEND example_options NO_TEST INSTALL_FILES idaRoberts_sps.out)
else()
  list(APPEND example_options EXAMPLE_TYPE develop)
endif()
sundials_add_example(
  idaRoberts_sps idaRoberts_sps.c
  ${example_options}
  LINK_LIBRARIES PRIVATE ${IDA_C_LIBS} sundials_nvecserial
                 sundials_sunlinsolsuperlumt ${SUPERLUMT_LIBRARIES}
                 ${EXE_EXTRA_LINK_LIBS})
