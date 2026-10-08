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
if(SUNDIALS_ENABLE_CVODE)
  set(_answer_files)
  foreach(_args
          "--np 2 2 --imex"
          "--np 2 2 --mri-arkstep"
          "--np 2 2 --mri-cvode-global"
          "--np 2 2 --mri-cvode-local")
    string(REPLACE " " "_" _suffix "${_args}")
    list(APPEND _answer_files
         ark_diffusion_reaction_p_${_suffix}.out)
  endforeach()
  foreach(_args
          "--np 2 2 --imex"
          "--np 2 2 --mri-arkstep"
          "--np 2 2 --mri-cvode-global"
          "--np 2 2 --mri-cvode-local")
    string(REPLACE " " "_" _suffix "${_args}")
    set(_target ark_diffusion_reaction_p)
    if(NOT TARGET ${_target})
      sundials_add_example(${_target} ark_diffusion_reaction_p.cpp NO_TEST
                            INSTALL_FILES ${_answer_files}
                            INCLUDE_DIRECTORIES PRIVATE ${PROJECT_SOURCE_DIR}/examples/utilities
                            LINK_LIBRARIES PRIVATE MPI::MPI_CXX sundials_arkode
                            sundials_cvode sundials_nvecparallel
                            sundials_nvecmpiplusx ${EXE_EXTRA_LINK_LIBS})
    endif()
    sundials_add_test(${_target}_${_suffix} ${_target} TEST_ARGS ${_args}
                      MPI_NPROCS 4 ANSWER_DIR ${CMAKE_CURRENT_SOURCE_DIR}
                      ANSWER_FILE ${_target}_${_suffix}.out
                      EXAMPLE_TYPE exclude-single)
  endforeach()
endif()
