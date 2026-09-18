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
sundials_add_executable(cvAdvDiff_bnd_omp cvAdvDiff_bnd_omp.c)
set_target_properties(cvAdvDiff_bnd_omp PROPERTIES FOLDER "Examples")
target_link_libraries(cvAdvDiff_bnd_omp PRIVATE ${CVODE_C_LIBS} sundials_nvecopenmp
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvAdvDiff_bnd_omp cvAdvDiff_bnd_omp TEST_ARGS 4
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvAdvDiff_bnd_omp.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvAdvDiff_bnd_omp.c cvAdvDiff_bnd_omp.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cvAdvDiff/c-openmp)
endif()
