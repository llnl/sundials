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
sundials_add_example(ark_brusselator1D_task_local_nls_f2003
                      ark_brusselator1D_task_local_nls_f2003.f90
                      TEST_NAME ark_brusselator1D_task_local_nls_f2003_--monitor
                      TEST_ARGS --monitor MPI_NPROCS 4 FLOAT_PRECISION 2
                      ANSWER_FILE ark_brusselator1D_task_local_nls_f2003_--monitor.out
                      EXAMPLE_TYPE develop
                      INSTALL_FILES ark_brusselator1D_task_local_nls_f2003_--monitor_--explicit_--tf_3.out
                                     ark_brusselator1D_task_local_nls_f2003_--monitor_--global-nls.out
                      LINK_LIBRARIES PRIVATE MPI::MPI_Fortran
                      sundials_arkode sundials_farkode_mod
                      sundials_nvecparallel sundials_fnvecparallel_mod
                      sundials_nvecmpiplusx sundials_fnvecmpiplusx_mod
                      sundials_nvecmpimanyvector sundials_fnvecmpimanyvector_mod
                      ${EXE_EXTRA_LINK_LIBS})
foreach(_args
        "--monitor --explicit --tf 3"
        "--monitor --global-nls")
  string(REPLACE " " "_" _suffix "${_args}")
  sundials_add_test(ark_brusselator1D_task_local_nls_f2003_${_suffix}
                    ark_brusselator1D_task_local_nls_f2003 TEST_ARGS ${_args}
                    MPI_NPROCS 4 ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                    ANSWER_FILE ark_brusselator1D_task_local_nls_f2003_${_suffix}.out
                    EXAMPLE_TYPE develop)
endforeach()
