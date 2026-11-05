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
sundials_add_example(
  ark_brusselator_1d_imexmri
  SOURCES ark_brusselator_1d_imexmri.c
  TEST_NAME ark_brusselator_1d_imexmri_0_0.001
  TEST_ARGS 0 0.001
  ANSWER_FILE ark_brusselator_1d_imexmri_0_0.001.out
  EXAMPLE_TYPE exclude-single
  INSTALL_FILES
    plot_brusselator_1d.py
    ark_brusselator_1d_imexmri_2_0.001.out
    ark_brusselator_1d_imexmri_3_0.001.out
    ark_brusselator_1d_imexmri_4_0.001.out
    ark_brusselator_1d_imexmri_5_0.001.out
    ark_brusselator_1d_imexmri_6_0.001.out
    ark_brusselator_1d_imexmri_7_0.001.out
  LINK_LIBRARIES PRIVATE sundials_arkode sundials_nvecserial
                 sundials_sundomeigestpower ${EXE_EXTRA_LINK_LIBS})
foreach(_case 2 3 4 5 6 7)
  sundials_add_test(
    ark_brusselator_1d_imexmri_${_case}_0.001 ark_brusselator_1d_imexmri
    TEST_ARGS ${_case} 0.001
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ark_brusselator_1d_imexmri_${_case}_0.001.out
    EXAMPLE_TYPE exclude-single)
endforeach()
