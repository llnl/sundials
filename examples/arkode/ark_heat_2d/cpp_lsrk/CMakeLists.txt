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
# Preserve develop's single-precision test exclusion for these cases.
sundials_add_example(
  ark_heat_2d_lsrk
  SOURCES ark_heat_2d_lsrk.cpp
  EXAMPLE_TYPE exclude-single
  INSTALL_FILES plot_heat_2d.py
  LINK_LIBRARIES PRIVATE sundials_arkode sundials_nvecserial
                 ${EXE_EXTRA_LINK_LIBS})
