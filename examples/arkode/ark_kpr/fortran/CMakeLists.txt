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
sundials_add_example(ark_kpr_nestedmri_f2003 ark_kpr_nestedmri_f2003.f90
                      EXAMPLE_TYPE develop
                      LINK_LIBRARIES PRIVATE sundials_farkode_mod sundials_arkode
                      ${EXE_EXTRA_LINK_LIBS})

if(SUNDIALS_INDEX_SIZE MATCHES "64")
  if(NOT TARGET ark_kpr_mri_f2003)
    file(GLOB _answer_files
         "${CMAKE_CURRENT_SOURCE_DIR}/ark_kpr_mri_f2003*.out")
    sundials_add_example(ark_kpr_mri_f2003 ark_kpr_mri_f2003.f90 NO_TEST
                          INSTALL_FILES ${_answer_files}
                          LINK_LIBRARIES PRIVATE sundials_farkode_mod sundials_arkode
                          ${EXE_EXTRA_LINK_LIBS})
  endif()
  set(_args "" "0 0.002" "1 0.002" "2 0.005" "3 0.01" "4 0.002"
            "5 0.002" "6 0.005" "7 0.001" "8 0.001" "9 0.001")
  foreach(_arg ${_args})
    if(_arg STREQUAL "")
      set(_suffix ark_kpr_mri_f2003)
      set(_args_list)
    else()
      string(REPLACE " " "_" _suffix ark_kpr_mri_f2003_${_arg})
      set(_args_list ${_arg})
    endif()
    sundials_add_test(${_suffix} ark_kpr_mri_f2003 TEST_ARGS ${_args_list}
                      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                      ANSWER_FILE ${_suffix}.out EXAMPLE_TYPE develop)
  endforeach()
endif()
