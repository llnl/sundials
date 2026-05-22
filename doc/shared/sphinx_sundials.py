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
from sphinx.domains.changeset import VersionChange

class SundialsVersionChange(VersionChange):
    def run(self):
        # Join all arguments into the version string so the full text appears before the
        # colon, e.g., "7.3.0 (ARKODE 6.3.0)". This requires the description to be on a
        # different line from the version text.
        self.arguments[0] = " ".join(self.arguments)
        self.arguments = [self.arguments[0]]
        return super().run()

def setup(app: Sphinx):
    # Create new object type for CMake options
    app.add_object_type("cmakeoption", "cmakeop", "single: CMake options; %s")
    # Create new configuration value set in conf.py
    app.add_config_value("package_name", "", "env", types=[str])
    # Override the built-in directives
    app.add_directive("versionchanged", SundialsVersionChange, override=True)
    app.add_directive("versionadded", SundialsVersionChange, override=True)
    app.add_directive("deprecated", SundialsVersionChange, override=True)
