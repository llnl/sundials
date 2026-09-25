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
sundials_add_executable(cv_pendulum_dns cv_pendulum_dns.c)
set_target_properties(cv_pendulum_dns PROPERTIES FOLDER "Examples")
target_link_libraries(cv_pendulum_dns PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cv_pendulum_dns cv_pendulum_dns ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cv_pendulum_dns.out EXAMPLE_TYPE exclude-single)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cv_pendulum_dns.c cv_pendulum_dns.out plot_cv_pendulum.py
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cv_pendulum/c)
endif()
