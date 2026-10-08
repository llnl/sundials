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
if(NOT TARGET ark_lotka_volterra_ASA)
  sundials_add_example(ark_lotka_volterra_ASA ark_lotka_volterra_ASA.c NO_TEST
                        INSTALL_FILES ark_lotka_volterra_ASA_--check-freq_1.out
                                       ark_lotka_volterra_ASA_--check-freq_5.out
                        LINK_LIBRARIES PRIVATE sundials_arkode sundials_nvecserial
                        ${EXE_EXTRA_LINK_LIBS})
endif()
foreach(_freq 1 5)
  sundials_add_test(ark_lotka_volterra_ASA_--check-freq_${_freq}
                    ark_lotka_volterra_ASA TEST_ARGS --check-freq ${_freq}
                    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                    ANSWER_FILE ark_lotka_volterra_ASA_--check-freq_${_freq}.out
                    EXAMPLE_TYPE develop)
endforeach()
