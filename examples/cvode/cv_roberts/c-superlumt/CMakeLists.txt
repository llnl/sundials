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
sundials_add_executable(cv_roberts_sps cv_roberts_sps.c)
set_target_properties(cv_roberts_sps PROPERTIES FOLDER "Examples")
target_link_libraries(cv_roberts_sps PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      sundials_sunlinsolsuperlumt ${SUPERLUMT_LIBRARIES}
                      ${EXE_EXTRA_LINK_LIBS})
if(SUNDIALS_INDEX_SIZE MATCHES "32")
  sundials_add_test(cv_roberts_sps cv_roberts_sps ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                    ANSWER_FILE cv_roberts_sps.out EXAMPLE_TYPE develop)
endif()
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cv_roberts_sps.c cv_roberts_sps.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cv_roberts/c-superlumt)
endif()
