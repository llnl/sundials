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
if(SUNDIALS_INDEX_SIZE MATCHES "64")
  sundials_add_example(
    ida_heat_2d_kry_omp_f2003 ida_heat_2d_kry_omp_f2003.f90
    NO_TEST
    INSTALL_FILES ida_heat_2d_kry_omp_f2003_4.out
                  ida_heat_2d_kry_omp_f2003_8.out
    LINK_LIBRARIES PRIVATE ${IDA_C_LIBS} sundials_fida_mod
                   sundials_nvecopenmp sundials_fnvecopenmp_mod
                   ${EXE_EXTRA_LINK_LIBS}
    PROPERTIES Fortran_MODULE_DIRECTORY
               ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/ida_heat_2d_kry_omp_f2003.dir)
  foreach(num_threads 4 8)
    sundials_add_test(
      ida_heat_2d_kry_omp_f2003_${num_threads} ida_heat_2d_kry_omp_f2003
      TEST_ARGS ${num_threads}
      ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
      ANSWER_FILE ida_heat_2d_kry_omp_f2003_${num_threads}.out
      EXAMPLE_TYPE exclude)
  endforeach()
endif()
