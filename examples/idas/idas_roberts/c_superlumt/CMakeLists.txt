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
foreach(example idas_roberts_asai_sps idas_roberts_fsa_sps idas_roberts_sps)
  set(example_options EXAMPLE_TYPE exclude-single)
  if(example STREQUAL "idas_roberts_fsa_sps")
    set(answer_file idas_roberts_fsa_sps_-sensi_stg_t.out)
    list(
      APPEND
      example_options
      TEST_NAME
      idas_roberts_fsa_sps_-sensi_stg_t
      TEST_ARGS
      -sensi
      stg
      t)
  else()
    set(answer_file ${example}.out)
  endif()
  if(NOT SUNDIALS_INDEX_SIZE MATCHES "32")
    list(APPEND example_options NO_TEST INSTALL_FILES ${answer_file})
  endif()
  sundials_add_example(${example} SOURCES ${example}.c ${example_options})
  target_link_libraries(
    ${example}
    PRIVATE ${IDAS_C_LIBS} sundials_nvecserial sundials_sunlinsolsuperlumt
            ${SUPERLUMT_LIBRARIES} ${EXE_EXTRA_LINK_LIBS})
endforeach()
