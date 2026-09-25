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
sundials_add_executable(cv_adv_diff_bnd_l cv_adv_diff_bnd_l.c)
set_target_properties(cv_adv_diff_bnd_l PROPERTIES FOLDER "Examples")
target_link_libraries(cv_adv_diff_bnd_l PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                      sundials_sunlinsollapackband sundials_sunlinsollapackdense
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cv_adv_diff_bnd_l cv_adv_diff_bnd_l ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cv_adv_diff_bnd_l.out EXAMPLE_TYPE exclude-single)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cv_adv_diff_bnd_l.c cv_adv_diff_bnd_l.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cv_adv_diff/c-lapack)
endif()
