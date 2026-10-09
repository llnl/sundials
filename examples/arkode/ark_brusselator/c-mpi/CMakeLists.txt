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
sundials_add_example(ark_brusselator1D_task_local_nls
                      ark_brusselator1D_task_local_nls.c
                      TEST_NAME ark_brusselator1D_task_local_nls_--monitor
                      TEST_ARGS --monitor MPI_NPROCS 4
                      ANSWER_FILE ark_brusselator1D_task_local_nls_--monitor.out
                      EXAMPLE_TYPE develop
                      INSTALL_FILES plot_brusselator1D.py
                                     ark_brusselator1D_task_local_nls_--monitor_--explicit_--tf_3.out
                                     ark_brusselator1D_task_local_nls_--monitor_--global-nls.out
                      LINK_LIBRARIES PRIVATE sundials_arkode sundials_nvecparallel
                      sundials_nvecmpiplusx ${EXE_EXTRA_LINK_LIBS})
foreach(_args
        "--monitor --explicit --tf 3"
        "--monitor --global-nls")
  string(REPLACE " " "_" _suffix "${_args}")
  set(_type develop)
  sundials_add_test(ark_brusselator1D_task_local_nls_${_suffix}
                    ark_brusselator1D_task_local_nls TEST_ARGS ${_args}
                    MPI_NPROCS 4 ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                    ANSWER_FILE ark_brusselator1D_task_local_nls_${_suffix}.out
                    EXAMPLE_TYPE ${_type})
endforeach()
