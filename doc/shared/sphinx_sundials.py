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


def source_replacements_handler(app, docname, source):
    for old, new in app.config.source_replacements.items():
        source[0] = source[0].replace(old, new)


def setup(app: Sphinx):
    # Create new object type for CMake options
    app.add_object_type("cmakeoption", "cmakeop", "single: CMake options; %s")
    # Create new configuration value sets in conf.py
    app.add_config_value("package_name", "", "env", types=[str])
    app.add_config_value("source_replacements", {}, "env", types=[dict])
    # Replace key strings in source_replacements dictionary with corresponding values.
    # This allows for replacing strings anywhere e.g., inside ..code-block instead of
    # just parsed-literal using rst_epilog
    app.connect("source-read", source_replacements_handler)
