# ------------------------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------

# test source files
set(target "ark_dae_fixed_m")
sundials_add_executable(${target} ${target}.cpp)
set_target_properties(${target} PROPERTIES FOLDER "${target}")
target_include_directories(${target} PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/test)
target_link_libraries(${target} sundials_arkode sundials_nvecserial ${EXE_EXTRA_LINK_LIBS})

set(target "ida_dae_fixed_m")
sundials_add_executable(${target} ${target}.cpp)
set_target_properties(${target} PROPERTIES FOLDER "${target}")
target_include_directories(${target} PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/test)
target_link_libraries(${target} sundials_ida sundials_nvecserial ${EXE_EXTRA_LINK_LIBS})

# Create symlink to plot.py in the build directory
file(CREATE_LINK "${CMAKE_CURRENT_SOURCE_DIR}/plot_dae_fixed_m.py" "${CMAKE_CURRENT_BINARY_DIR}/plot_dae_fixed_m.py" SYMBOLIC)
