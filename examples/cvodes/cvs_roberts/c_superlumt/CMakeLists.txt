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
foreach(example cvs_roberts_asai_sps cvs_roberts_fsa_sps cvs_roberts_sps)
  set(answer_file ${example}.out)
  set(example_options)
  if(example STREQUAL "cvs_roberts_asai_sps")
    list(APPEND example_options EXAMPLE_TYPE exclude-single)
  elseif(example STREQUAL "cvs_roberts_fsa_sps")
    set(answer_file cvs_roberts_fsa_sps_-sensi_stg1_t.out)
    list(
      APPEND
      example_options
      TEST_NAME
      cvs_roberts_fsa_sps_-sensi_stg1_t
      TEST_ARGS
      -sensi
      stg1
      t
      EXAMPLE_TYPE
      exclude-single)
  endif()
  if(NOT SUNDIALS_INDEX_SIZE MATCHES "32")
    list(APPEND example_options NO_TEST INSTALL_FILES ${answer_file})
  endif()
  sundials_add_example(${example} ${example}.c ${example_options})
  target_link_libraries(
    ${example}
    PRIVATE ${CVODES_C_LIBS} sundials_nvecserial sundials_sunlinsolsuperlumt
            ${SUPERLUMT_LIBRARIES} ${EXE_EXTRA_LINK_LIBS})
endforeach()
