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
set(_xbraid_examples ark_heat2D_xbraid ark_heat2D_p_xbraid)
if(SUNDIALS_ENABLE_HYPRE)
  list(APPEND _xbraid_examples ark_heat2D_hypre_pfmg_xbraid)
endif()
foreach(_example ${_xbraid_examples})
  if(_example STREQUAL "ark_heat2D_xbraid")
    set(_args "--x_print_level 0")
    set(_tasks 4)
    set(_type develop)
  elseif(_example STREQUAL "ark_heat2D_p_xbraid")
    set(_args "--np 2 1 2 --x_print_level 0")
    set(_tasks 4)
    set(_type develop)
  else()
    set(_args "--np 2 1 2 --x_print_level 0")
    set(_tasks 4)
    set(_type develop)
  endif()
  string(REPLACE " " "_" _suffix "${_args}")
  set(_libs MPI::MPI_CXX sundials_arkode sundials_arkode_xbraid
            sundials_nvecparallel ${EXE_EXTRA_LINK_LIBS})
  if(_example MATCHES "hypre")
    list(APPEND _libs SUNDIALS::HYPRE SUNDIALS::XBRAID)
  endif()
  sundials_add_example(${_example} ${_example}.cpp NO_TEST
                        INSTALL_FILES plot_heat2D.py ${_example}_${_suffix}.out
                        LINK_LIBRARIES PRIVATE ${_libs})
  sundials_add_test(${_example}_${_suffix} ${_example} TEST_ARGS ${_args}
                    MPI_NPROCS ${_tasks} ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                    ANSWER_FILE ${_example}_${_suffix}.out
                    EXAMPLE_TYPE ${_type})
endforeach()
