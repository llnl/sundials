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
sundials_add_executable(cv_roberts_dns cv_roberts_dns.c)
set_target_properties(cv_roberts_dns PROPERTIES FOLDER "Examples")
target_link_libraries(cv_roberts_dns PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cv_roberts_dns cv_roberts_dns ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cv_roberts_dns.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cv_roberts_dns.c cv_roberts_dns.out cv_roberts_dns_stats.csv
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cvRoberts/c)
endif()

foreach(example cv_roberts_dns_constraints cv_roberts_dns_negsol cv_roberts_dns_uw)
  sundials_add_executable(${example} ${example}.c)
  set_target_properties(${example} PROPERTIES FOLDER "Examples")
  target_link_libraries(${example} PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                        ${EXE_EXTRA_LINK_LIBS})
  sundials_add_test(${example} ${example} ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                    ANSWER_FILE ${example}.out EXAMPLE_TYPE develop)
endforeach()
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cv_roberts_dns_constraints.c cv_roberts_dns_constraints.out
                cv_roberts_dns_negsol.c cv_roberts_dns_negsol.out
                cv_roberts_dns_uw.c cv_roberts_dns_uw.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cvRoberts/c)
endif()
