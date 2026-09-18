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
sundials_add_executable(cvRoberts_dnsL cvRoberts_dnsL.c)
set_target_properties(cvRoberts_dnsL PROPERTIES FOLDER "Examples")
target_link_libraries(cvRoberts_dnsL PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      sundials_sunlinsollapackband sundials_sunlinsollapackdense
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvRoberts_dnsL cvRoberts_dnsL ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvRoberts_dnsL.out EXAMPLE_TYPE exclude-single)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvRoberts_dnsL.c cvRoberts_dnsL.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cvRoberts/c-lapack)
endif()
