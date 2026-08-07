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

# Steps to add an external addon to SUNDIALS:
#
# WARNING: external addons are not maintained by the SUNDIALS team and may come
# with their own license with different rules.
#
# 1. Clone/copy the addon into <sundials root>/external/ (the same directory as
#    this file you are reading).
# 2. Copy the sundials-addon-example block for your own external project.
# 3. When building SUNDIALS, set the CMake option SUNDIALS_ENABLE_EXTERNAL_ADDONS
#    to ON
# 4. Build SUNDIALS as usual.

include(FetchContent)

# COPY THE FetchContent BLOCK BELOW FOR YOUR OWN EXTERNAL ADDON COMMENT THESE
# LINES OUT TO DISABLE THE EXAMPLE ADDON
FetchContent_Declare(
  sundials-addon-example SOURCE_DIR
                         ${PROJECT_SOURCE_DIR}/external/sundials-addon-example)
FetchContent_MakeAvailable(sundials-addon-example)
# COPY THE BLOCK ABOVE FOR YOUR OWN EXTERNAL ADDON
