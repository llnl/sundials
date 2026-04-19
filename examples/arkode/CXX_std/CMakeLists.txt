# ---------------------------------------------------------------
# Programmer(s): David J. Gardner @ LLNL
# ---------------------------------------------------------------
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
# ---------------------------------------------------------------

add_executable(harmonic harmonic_oscillator.cpp)
target_link_libraries(harmonic PRIVATE
  sundials_arkode sundials_nvecserial
)

add_executable(harmonic_serial harmonic_oscillator_serial.cpp)
target_link_libraries(harmonic_serial PRIVATE
  sundials_arkode sundials_nvecserial
)

add_executable(coupled_oscillators coupled_oscillators.cpp)
target_link_libraries(coupled_oscillators PRIVATE
  sundials_arkode sundials_nvecserial
)

add_executable(coupled_oscillators_serial coupled_oscillators_serial.cpp)
target_link_libraries(coupled_oscillators_serial PRIVATE
  sundials_arkode sundials_nvecserial
)
