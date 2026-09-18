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
# CVODE parallel C++ heat equation example
set(CVODE_examples "cv_heat2D_p.cpp\;--np 2 2\;1\;4\;develop")
sundials_add_executable(cv_heat2D_p cv_heat2D_p.cpp)
set_target_properties(cv_heat2D_p PROPERTIES FOLDER "Examples")
target_link_libraries(cv_heat2D_p PRIVATE MPI::MPI_CXX sundials_cvode
                      sundials_nvecparallel ${EXE_EXTRA_LINK_LIBS})
if(SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS)
  target_link_libraries(cv_heat2D_p PRIVATE sundials_cvode_fused_stubs)
endif()
sundials_add_test(cv_heat2D_p_--np_2_2 cv_heat2D_p TEST_ARGS --np 2 2 MPI_NPROCS 4
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cv_heat2D_p_--np_2_2.out EXAMPLE_TYPE develop)

if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cv_heat2D_p.cpp cv_heat2D_p_--np_2_2.out plot_heat2D_p.py
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cvHeat2D/cpp-mpi)
endif()
