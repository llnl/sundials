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
sundials_add_executable(cvAdvDiff_non_ph cvAdvDiff_non_ph.c)
set_target_properties(cvAdvDiff_non_ph PROPERTIES FOLDER "Examples")
target_link_libraries(cvAdvDiff_non_ph PRIVATE ${CVODE_C_LIBS} sundials_nvecparhyp
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cvAdvDiff_non_ph cvAdvDiff_non_ph MPI_NPROCS 2
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cvAdvDiff_non_ph.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cvAdvDiff_non_ph.c cvAdvDiff_non_ph.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cvAdvDiff/c-mpi-hypre)
endif()
