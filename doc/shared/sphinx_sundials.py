# -----------------------------------------------------------------------------
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
# -----------------------------------------------------------------------------
# SUNDIALS Sphinx extension
# -----------------------------------------------------------------------------

from sphinx.application import Sphinx


def setup(app: Sphinx):
    # Create new object type for CMake options
    app.add_object_type("cmakeoption", "cmakeop", "single: CMake options; %s")
    # Create new configuration value set in conf.py
    app.add_config_value("package_name", "", "env", types=[str])
