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
# Tests are registered separately below for their argument-specific answers.
set(_args "--x_print_level 0")
set(_tasks 4)
string(REPLACE " " "_" _suffix "${_args}")
set(_libs MPI::MPI_CXX sundials_arkode sundials_arkode_xbraid
          sundials_nvecparallel ${EXE_EXTRA_LINK_LIBS})
sundials_add_example(
  ark_heat_2d_xbraid
  SOURCES ark_heat_2d_xbraid.cpp
  NO_TEST
  INSTALL_FILES plot_heat_2d.py ark_heat_2d_xbraid_${_suffix}.out
  LINK_LIBRARIES PRIVATE ${_libs})
sundials_add_test(
  ark_heat_2d_xbraid_${_suffix} ark_heat_2d_xbraid
  TEST_ARGS ${_args}
  MPI_NPROCS ${_tasks}
  ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  ANSWER_FILE ark_heat_2d_xbraid_${_suffix}.out
  EXAMPLE_TYPE develop)
