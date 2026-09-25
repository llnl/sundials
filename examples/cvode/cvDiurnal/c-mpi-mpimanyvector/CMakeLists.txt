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
sundials_add_executable(cv_diurnal_kry_mpimanyvec cv_diurnal_kry_mpimanyvec.c)
set_target_properties(cv_diurnal_kry_mpimanyvec PROPERTIES FOLDER "Examples")
target_link_libraries(cv_diurnal_kry_mpimanyvec PRIVATE MPI::MPI_C ${CVODE_C_LIBS}
                      sundials_nvecmpimanyvector sundials_nvecparallel
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_test(cv_diurnal_kry_mpimanyvec cv_diurnal_kry_mpimanyvec MPI_NPROCS 4
                  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE cv_diurnal_kry_mpimanyvec.out EXAMPLE_TYPE develop)
if(SUNDIALS_ENABLE_EXAMPLES_INSTALL)
  install(FILES cv_diurnal_kry_mpimanyvec.c cv_diurnal_kry_mpimanyvec.out
          DESTINATION ${SUNDIALS_EXAMPLES_INSTALL_PATH}/cvode/cvDiurnal/c-mpi-mpimanyvector)
endif()
