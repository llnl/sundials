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
foreach(_example ark_heat2D_p ark_heat2D_lsrk_p)
  set(_type develop)
  if(_example STREQUAL "ark_heat2D_lsrk_p")
    set(_type exclude-single)
  endif()
  set(_args "--np 2 2")
  string(REPLACE " " "_" _suffix "${_args}")
  sundials_add_example(${_example} ${_example}.cpp TEST_ARGS ${_args}
                        MPI_NPROCS 4
                        ANSWER_FILE ${_example}_${_suffix}.out
                        EXAMPLE_TYPE ${_type}
                        INSTALL_FILES plot_heat2D_p.py plot_data_2d.py
                        LINK_LIBRARIES PRIVATE MPI::MPI_CXX sundials_arkode
                        sundials_nvecparallel ${EXE_EXTRA_LINK_LIBS})
endforeach()
