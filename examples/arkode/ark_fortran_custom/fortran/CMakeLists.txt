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
set(_sources fnvector_complex_mod.f90 fnvector_fortran_mod.f90
             fsunmatrix_fortran_mod.f90 fsunlinsol_fortran_mod.f90)

sundials_add_example(ark_analytic_complex_f2003 ark_analytic_complex_f2003.f90
                      ${_sources} EXAMPLE_TYPE develop
                      INSTALL_FILES ${_sources}
                      PROPERTIES Fortran_MODULE_DIRECTORY
                                 ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/ark_analytic_complex_f2003.dir
                      LINK_LIBRARIES PRIVATE sundials_farkode_mod sundials_arkode
                      ${EXE_EXTRA_LINK_LIBS})
sundials_add_example(ark_brusselator1D_f2003 ark_brusselator1D_f2003.f90
                      ${_sources} EXAMPLE_TYPE develop
                      INSTALL_FILES ${_sources}
                      PROPERTIES Fortran_MODULE_DIRECTORY
                                 ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/ark_brusselator1D_f2003.dir
                      LINK_LIBRARIES PRIVATE sundials_farkode_mod sundials_arkode
                      ${EXE_EXTRA_LINK_LIBS})

foreach(_test test_fnvector_complex_mod test_fnvector_fortran_mod
              test_fsunmatrix_fortran_mod test_fsunlinsol_fortran_mod)
  sundials_add_example(${_test} ${_test}.f90 ${_sources} NODIFF
                        EXAMPLE_TYPE develop INSTALL_FILES ${_sources}
                        PROPERTIES Fortran_MODULE_DIRECTORY
                                   ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${_test}.dir
                        LINK_LIBRARIES PRIVATE sundials_farkode_mod sundials_arkode
                        ${EXE_EXTRA_LINK_LIBS})
endforeach()
