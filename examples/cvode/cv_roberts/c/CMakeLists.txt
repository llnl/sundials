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
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  cmake_minimum_required(VERSION 3.18)
  project(cvode_example LANGUAGES NONE)
endif()
include(${CMAKE_CURRENT_LIST_DIR}/../../SundialsCVODEExample.cmake)
foreach(example cv_roberts_dns cv_roberts_dns_constraints cv_roberts_dns_negsol
                cv_roberts_dns_uw)
  sundials_add_executable(${example} ${example}.c)
  set_target_properties(${example} PROPERTIES FOLDER "Examples")
  target_link_libraries(${example} PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                                           ${EXE_EXTRA_LINK_LIBS})
  if(example STREQUAL "cv_roberts_dns_negsol")
    set(example_type exclude-single)
  else()
    set(example_type develop)
  endif()
  sundials_add_test(
    ${example} ${example}
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}.out
    EXAMPLE_TYPE ${example_type})
endforeach()

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(
    FILES cv_roberts_dns.c
          cv_roberts_dns.out
          cv_roberts_dns_stats.csv
          cv_roberts_dns_constraints.c
          cv_roberts_dns_constraints.out
          cv_roberts_dns_negsol.c
          cv_roberts_dns_negsol.out
          cv_roberts_dns_uw.c
          cv_roberts_dns_uw.out
    DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cv_roberts/c)
endif()
