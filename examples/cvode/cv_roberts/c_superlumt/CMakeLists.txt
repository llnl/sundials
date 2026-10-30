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
if(SUNDIALS_INDEX_SIZE MATCHES "32")
  set(_no_test "")
else()
  set(_no_test NO_TEST)
endif()

sundials_add_example(cv_roberts_sps cv_roberts_sps.c
  ${_no_test}
  LINK_LIBRARIES
  PRIVATE ${CVODE_C_LIBS} sundials_nvecserial sundials_sunlinsolsuperlumt
          ${SUPERLUMT_LIBRARIES} ${EXE_EXTRA_LINK_LIBS})
