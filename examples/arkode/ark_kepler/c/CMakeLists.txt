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
    "--stepper ERK --step-mode adapt"
    "--stepper ERK --step-mode fixed --count-orbits"
    "--stepper SPRK --step-mode fixed --count-orbits --use-compensated-sums"
    "--stepper SPRK --step-mode fixed --method ARKODE_SPRK_EULER_1_1 --tf 50 --check-order --nout 1"
    "--stepper SPRK --step-mode fixed --method ARKODE_SPRK_LEAPFROG_2_2 --tf 50 --check-order --nout 1"
    "--stepper SPRK --step-mode fixed --method ARKODE_SPRK_MCLACHLAN_2_2 --tf 50 --check-order --nout 1"
    "--stepper SPRK --step-mode fixed --method ARKODE_SPRK_MCLACHLAN_3_3 --tf 50 --check-order --nout 1"
    "--stepper SPRK --step-mode fixed --method ARKODE_SPRK_MCLACHLAN_4_4 --tf 50 --check-order --nout 1"
    "--stepper SPRK --step-mode fixed --method ARKODE_SPRK_MCLACHLAN_5_6 --tf 50 --check-order --nout 1"
    "--stepper SPRK --step-mode fixed --method ARKODE_SPRK_PSEUDO_LEAPFROG_2_2 --tf 50 --check-order --nout 1"
    "--stepper SPRK --step-mode fixed --method ARKODE_SPRK_RUTH_3_3 --tf 50 --check-order --nout 1"
    "--stepper SPRK --step-mode fixed --method ARKODE_SPRK_YOSHIDA_6_8 --tf 50 --check-order --nout 1")

set(_answer_files ark_kepler.out)
foreach(_arg ${_args})
  string(REPLACE " " "_" _suffix "${_arg}")
  list(APPEND _answer_files ark_kepler_${_suffix}.out)
endforeach()

if(NOT TARGET ark_kepler)
  sundials_add_example(ark_kepler ark_kepler.c NO_TEST
                        INSTALL_FILES ark_kepler.h ark_kepler_plot.py
                                       ${_answer_files}
                        LINK_LIBRARIES PRIVATE sundials_arkode
                        sundials_nvecserial ${EXE_EXTRA_LINK_LIBS})
endif()
sundials_add_test(ark_kepler ark_kepler ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                  ANSWER_FILE ark_kepler.out EXAMPLE_TYPE develop)

foreach(_arg ${_args})
  string(REPLACE " " "_" _suffix "${_arg}")
  set(_type develop)
  if(_arg MATCHES "check-order")
    set(_type exclude-single)
  endif()
  sundials_add_test(ark_kepler_${_suffix} ark_kepler TEST_ARGS ${_arg}
                    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                    ANSWER_FILE ark_kepler_${_suffix}.out
                    EXAMPLE_TYPE ${_type})
endforeach()
