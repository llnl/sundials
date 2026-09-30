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
# CVODE SuperLU_DIST C++ example
sundials_add_example(
  cv_adv_diff_sludist cv_adv_diff_sludist.cpp
  TEST_ARGS 2 1
  MPI_NPROCS 2
  ANSWER_FILE cv_adv_diff_sludist_${SUNDIALS_INDEX_SIZE}.out
  INTEGER_PRECISION ${_integer_precision}
  FLOAT_PRECISION ${_float_precision}
  INSTALL_FILES cv_adv_diff_sludist_32.out cv_adv_diff_sludist_64.out)
target_link_libraries(
  cv_adv_diff_sludist
  PRIVATE MPI::MPI_CXX ${CVODE_C_LIBS} sundials_nvecparallel
          sundials_sunlinsolsuperludist ${EXE_EXTRA_LINK_LIBS})
if(SUNDIALS_INDEX_SIZE STREQUAL "32")
  set(_integer_precision 15)
  set(_float_precision 3)
else()
  set(_integer_precision default)
  set(_float_precision default)
endif()
