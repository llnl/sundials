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
# CVODE serial C++ heat equation example
sundials_add_example(
  cv_heat_2d cv_heat_2d.cpp
  EXAMPLE_TYPE exclude-single
  INSTALL_FILES cv_heat_2d.hpp plot_heat_2d.py
                "${PROJECT_SOURCE_DIR}/examples/utilities/example_utilities.hpp"
  INCLUDE_DIRECTORIES PRIVATE ${PROJECT_SOURCE_DIR}/examples/utilities
  LINK_LIBRARIES PRIVATE ${CVODE_C_LIBS} sundials_nvecserial
                 ${EXE_EXTRA_LINK_LIBS})
