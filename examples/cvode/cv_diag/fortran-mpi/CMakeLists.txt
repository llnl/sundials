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
foreach(example cv_diag_kry_bbd_f2003 cv_diag_kry_f2003 cv_diag_non_p_f2003)
  set(CMAKE_Fortran_MODULE_DIRECTORY
      ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${example}.dir)
  sundials_add_example(
    ${example} ${example}.f90 ADD_TEST
    MPI_NPROCS 4
    ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    ANSWER_FILE ${example}.out
    EXAMPLE_TYPE develop
    FLOAT_PRECISION 2
    INSTALL_FILES ${example}.out)
  set_target_properties(${example} PROPERTIES FOLDER "Examples")
  target_link_libraries(
    ${example}
    MPI::MPI_Fortran
    sundials_cvode
    sundials_fcvode_mod
    sundials_nvecparallel
    sundials_fnvecparallel_mod
    ${EXE_EXTRA_LINK_LIBS})
endforeach()
