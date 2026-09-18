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
sundials_add_executable(cvRoberts_dns cvRoberts_dns.c)
set_target_properties(cvRoberts_dns PROPERTIES FOLDER "Examples")
target_link_libraries(cvRoberts_dns PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvRoberts_dns cvRoberts_dns ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvRoberts_dns.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvRoberts_dns.c cvRoberts_dns.out cvRoberts_dns_stats.csv
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cvRoberts/c)
endif()

foreach(example cvRoberts_dns_constraints cvRoberts_dns_negsol cvRoberts_dns_uw)
  sundials_add_executable(${example} ${example}.c)
  set_target_properties(${example} PROPERTIES FOLDER "Examples")
  target_link_libraries(${example} PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                        ${EXE_EXTRA_LINK_LIBS})
  sundials_add_test(${example} ${example} ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                    ANSWER_FILE ${example}.out EXAMPLE_TYPE develop)
endforeach()
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvRoberts_dns_constraints.c cvRoberts_dns_constraints.out
                cvRoberts_dns_negsol.c cvRoberts_dns_negsol.out
                cvRoberts_dns_uw.c cvRoberts_dns_uw.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cvRoberts/c)
endif()
