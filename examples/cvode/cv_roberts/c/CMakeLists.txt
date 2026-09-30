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
foreach(example cv_roberts_dns cv_roberts_dns_constraints cv_roberts_dns_negsol
                cv_roberts_dns_uw)
  if(example STREQUAL "cv_roberts_dns_negsol")
    set(example_type exclude-single)
  else()
    set(example_type develop)
  endif()
  set(install_files ${example}.out)
  if(example STREQUAL "cv_roberts_dns")
    list(APPEND install_files cv_roberts_dns_stats.csv)
  endif()
  sundials_add_example(
    ${example} ${example}.c ADD_TEST
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}.out
    EXAMPLE_TYPE ${example_type}
    INSTALL_FILES ${install_files})
  set_target_properties(${example} PROPERTIES FOLDER "Examples")
  target_link_libraries(${example} PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                                           ${EXE_EXTRA_LINK_LIBS})
endforeach()
