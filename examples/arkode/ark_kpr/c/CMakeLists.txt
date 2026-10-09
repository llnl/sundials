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
set(_args
    "0 1 0.005" "1 0 0.01" "1 1 0.002" "2 4 0.002" "3 2 0.001"
    "4 3 0.001" "5 4 0.001" "6 5 0.001" "7 2 0.002"
    "8 3 0.001 -100 100 0.5 1" "9 3 0.001 -100 100 0.5 1"
    "10 4 0.001 -100 100 0.5 1" "11 2 0.001" "12 3 0.005"
    "13 4 0.01")
file(GLOB _answer_files "${CMAKE_CURRENT_SOURCE_DIR}/ark_kpr_mri_*.out")
if(NOT TARGET ark_kpr_mri)
  sundials_add_example(ark_kpr_mri ark_kpr_mri.c NO_TEST
                        INSTALL_FILES ${_answer_files}
                        LINK_LIBRARIES PRIVATE sundials_arkode sundials_nvecserial
                        ${EXE_EXTRA_LINK_LIBS})
endif()
foreach(_arg ${_args})
  string(REPLACE " " "_" _suffix "${_arg}")
  set(_type develop)
  if(_arg MATCHES "^1 1|^8 |^9 |^10 |^13 ")
    set(_type exclude-single)
  endif()
  sundials_add_test(ark_kpr_mri_${_suffix} ark_kpr_mri TEST_ARGS ${_arg}
                    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                    ANSWER_FILE ark_kpr_mri_${_suffix}.out
                    EXAMPLE_TYPE ${_type})
endforeach()
