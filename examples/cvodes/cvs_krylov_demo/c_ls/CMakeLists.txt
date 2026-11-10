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
  cvs_krylov_demo_ls
  SOURCES cvs_krylov_demo_ls.c
  EXAMPLE_TYPE develop
  INSTALL_FILES cvs_krylov_demo_ls_0_1.out cvs_krylov_demo_ls_1.out
                cvs_krylov_demo_ls_2.out
  LINK_LIBRARIES ${CVODES_C_LIBS} sundials_nvecserial ${EXE_EXTRA_LINK_LIBS})
set(_test_args 1 2)
if(SUNDIALS_ENABLE_MONITORING)
  list(APPEND _test_args "0 1")
endif()
foreach(test_args ${_test_args})
  string(REPLACE " " "_" test_suffix "${test_args}")
  set(_answer_file cvs_krylov_demo_ls_${test_suffix}.out)
  if(test_args STREQUAL "0 1")
    # Monitoring writes to a log file; stdout matches the default run.
    set(_answer_file cvs_krylov_demo_ls.out)
  endif()
  sundials_add_test(
    cvs_krylov_demo_ls_${test_suffix} cvs_krylov_demo_ls
    TEST_ARGS ${test_args}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${_answer_file}
    EXAMPLE_TYPE develop)
endforeach()
