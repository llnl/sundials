# -----------------------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025, Lawrence Livermore National Security,
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
# Script that generates sphinx-autodoc autofunction directives for each
# function in the Python modules. This is necessary due to automodule not
# yet being able to handle nanobind generated functions.
# See https://github.com/sphinx-doc/sphinx/issues/13868.
# -----------------------------------------------------------------------------

import os
import importlib
from sphinx.util import logging

logger = logging.getLogger(__name__)

if os.getenv("SPHINX_MOCK_SUNDIALS4PY", "").lower() in ("true", "1", "yes", "on"):
    import sys
    from unittest.mock import MagicMock

    logger.info("Using mock sundials4py module.")

    class Mock(MagicMock):
        @classmethod
        def __getattr__(cls, name):
            if name == "__all__":
                return []  # Return empty list for __all__
            return MagicMock()

    # Mock sundials4py and submodules
    mock_module = Mock()
    mock_module.__all__ = []

    sys.modules["sundials4py"] = mock_module
    sys.modules["sundials4py.core"] = Mock()
    sys.modules["sundials4py.core"].__all__ = []
    sys.modules["sundials4py.kinsol"] = Mock()
    sys.modules["sundials4py.kinsol"].__all__ = []
    sys.modules["sundials4py.cvode"] = Mock()
    sys.modules["sundials4py.cvode"].__all__ = []
    sys.modules["sundials4py.cvodes"] = Mock()
    sys.modules["sundials4py.cvodes"].__all__ = []
    sys.modules["sundials4py.arkode"] = Mock()
    sys.modules["sundials4py.arkode"].__all__ = []
    sys.modules["sundials4py.ida"] = Mock()
    sys.modules["sundials4py.ida"].__all__ = []
    sys.modules["sundials4py.idas"] = Mock()
    sys.modules["sundials4py.idas"].__all__ = []

else:

    try:
        import sundials4py
    except ModuleNotFoundError as e:
        logger.error(
            f"{e}. Either install sundials4py or set the environment variable SPHINX_MOCK_SUNDIALS4PY to ON."
        )
        raise  # Re-raises the same exception


def generate_autofunctions_for_submodule(module_name: str):
    module = importlib.import_module(f"sundials4py.{module_name}")
    autogen_file = os.path.join(
        os.path.dirname(__file__), f"./Python/sundials4py-{module_name}-functions.rst"
    )
    with open(autogen_file, "w") as f:
        header = f"{module_name} Submodule"
        f.write(f"{header}\n")
        f.write(f"{('').join('=' for i in range(len(header)))}\n\n")
        f.write("Classes\n")
        f.write("^^^^^^^\n\n")
        f.write(f".. automodule:: sundials4py.{module_name}\n")
        f.write("   :members:\n")
        f.write("   :undoc-members:\n")
        f.write("   :private-members:\n\n")
        f.write("Functions\n")
        f.write("^^^^^^^^^\n\n")
        for func_name in dir(module):
            obj = getattr(module, func_name)
            if type(obj).__name__ == "nb_func":
                f.write(f".. autofunction:: sundials4py.{module_name}.{func_name}\n")
                f.write("  :no-index:\n\n")
                f.write(f"  See :c:func:`{func_name}`.\n\n")


def generate_autofunctions_for_sundials4py():
    generate_autofunctions_for_submodule("core")
    generate_autofunctions_for_submodule("arkode")
    generate_autofunctions_for_submodule("cvodes")
    generate_autofunctions_for_submodule("idas")
    generate_autofunctions_for_submodule("kinsol")
