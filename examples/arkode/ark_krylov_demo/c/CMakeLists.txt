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
sundials_add_example(ark_KrylovDemo_prec ark_KrylovDemo_prec.c
                      ANSWER_FILE ark_KrylovDemo_prec.out
                      EXAMPLE_TYPE exclude-single
                      INSTALL_FILES ark_KrylovDemo_prec_1.out
                                     ark_KrylovDemo_prec_2.out
                      LINK_LIBRARIES PRIVATE sundials_arkode sundials_nvecserial
                      ${EXE_EXTRA_LINK_LIBS})
foreach(_arg 1 2)
  sundials_add_test(ark_KrylovDemo_prec_${_arg} ark_KrylovDemo_prec
                    TEST_ARGS ${_arg} ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                    ANSWER_FILE ark_KrylovDemo_prec_${_arg}.out
                    EXAMPLE_TYPE exclude-single)
endforeach()
