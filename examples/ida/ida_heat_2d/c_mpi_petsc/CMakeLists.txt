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
  ida_heat_2d_petsc_spgmr ida_heat_2d_petsc_spgmr.c
  MPI_NPROCS 4
  EXAMPLE_TYPE develop
  LINK_LIBRARIES PRIVATE MPI::MPI_C ${IDA_C_LIBS} sundials_nvecpetsc
                 sundials_sunnonlinsolpetscsnes ${EXE_EXTRA_LINK_LIBS})

sundials_add_example(
  ida_heat_2d_petsc_snes ida_heat_2d_petsc_snes.c
  MPI_NPROCS 4
  EXAMPLE_TYPE develop
  INSTALL_FILES ida_heat_2d_petsc_snes_-pre.out ida_heat_2d_petsc_snes_-jac.out
                ida_heat_2d_petsc_snes_-jac_-pre.out
  LINK_LIBRARIES PRIVATE MPI::MPI_C ${IDA_C_LIBS} sundials_nvecpetsc
                 sundials_sunnonlinsolpetscsnes ${EXE_EXTRA_LINK_LIBS})

sundials_add_test(
  ida_heat_2d_petsc_snes_-pre ida_heat_2d_petsc_snes
  TEST_ARGS -pre
  MPI_NPROCS 4
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE ida_heat_2d_petsc_snes_-pre.out
  EXAMPLE_TYPE develop)
if(SUNDIALS_INDEX_SIZE STREQUAL "32")
  foreach(test_args "-jac" "-jac -pre")
    string(REPLACE " " "_" test_suffix "${test_args}")
    sundials_add_test(
      ida_heat_2d_petsc_snes_${test_suffix} ida_heat_2d_petsc_snes
      TEST_ARGS ${test_args}
      MPI_NPROCS 4
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ida_heat_2d_petsc_snes_${test_suffix}.out
      EXAMPLE_TYPE develop)
  endforeach()
endif()
