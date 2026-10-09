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
set(_arkode_libs sundials_arkode sundials_nvecserial
                 sundials_sundomeigestpower ${EXE_EXTRA_LINK_LIBS})

foreach(_example ark_brusselator ark_brusselator1D ark_brusselator1D_imexmri
                 ark_brusselator_1D_mri ark_brusselator_fp
                 ark_brusselator_lsrk_domeigest ark_brusselator_lsrk_externaldomeigest
                 ark_brusselator_mri)
  set(_type develop)
  if(_example MATCHES "1D$|imexmri|fp$")
    set(_type exclude-single)
  endif()
  if(_example STREQUAL "ark_brusselator_fp")
    sundials_add_example(${_example} ${_example}.c EXAMPLE_TYPE ${_type}
                          INSTALL_FILES ark_brusselator_fp_1.out
                          LINK_LIBRARIES PRIVATE ${_arkode_libs})
    if(SUNDIALS_ENABLE_MONITORING)
      sundials_add_test(ark_brusselator_fp_1 ark_brusselator_fp TEST_ARGS 1
                        ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                        ANSWER_FILE ark_brusselator_fp_1.out
                        EXAMPLE_TYPE exclude-single)
    endif()
  elseif(_example STREQUAL "ark_brusselator1D_imexmri")
    sundials_add_example(${_example} ${_example}.c TEST_NAME ${_example}_0_0.001
                          TEST_ARGS 0 0.001
                          ANSWER_FILE ${_example}_0_0.001.out
                          EXAMPLE_TYPE exclude-single
                          INSTALL_FILES ark_brusselator1D_imexmri_2_0.001.out
                                         ark_brusselator1D_imexmri_3_0.001.out
                                         ark_brusselator1D_imexmri_4_0.001.out
                                         ark_brusselator1D_imexmri_5_0.001.out
                                         ark_brusselator1D_imexmri_6_0.001.out
                                         ark_brusselator1D_imexmri_7_0.001.out
                          LINK_LIBRARIES PRIVATE ${_arkode_libs})
    foreach(_case 2 3 4 5 6 7)
      sundials_add_test(${_example}_${_case}_0.001 ${_example}
                        TEST_ARGS ${_case} 0.001
                        ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                        ANSWER_FILE ${_example}_${_case}_0.001.out
                        EXAMPLE_TYPE exclude-single)
    endforeach()
  else()
    sundials_add_example(${_example} ${_example}.c EXAMPLE_TYPE ${_type}
                          LINK_LIBRARIES PRIVATE ${_arkode_libs})
  endif()
endforeach()
