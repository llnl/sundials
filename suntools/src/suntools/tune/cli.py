#!/usr/bin/env python3
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

from __future__ import annotations

import sys
from typing import Any

from pydantic import ValidationError

from suntools.tune.config import config_from_args
from suntools.tune.deephyper_backend import DeepHyperBackend
from suntools.tune.gptune_backend import GPTuneBackend
from suntools.tune.runner import format_command, select_best
from suntools.tune.ytopt_backend import YtoptBackend


def run_from_args(args: Any) -> int:
    try:
        config = config_from_args(args)
        backend = _create_backend(config)
        results = backend.run()
    except (RuntimeError, ValueError, ValidationError) as err:
        sys.stderr.write("error: %s\n" % err)
        return 2

    best = select_best(results, config.objective.direction)
    if best is None:
        sys.stderr.write("error: no successful tune trials completed\n")
        return 1

    sys.stdout.write("Best %s: %s\n" % (config.objective.metric, best.metric))
    for name in sorted(best.parameters):
        sys.stdout.write("  %s = %s\n" % (name, best.parameters[name]))
    sys.stdout.write("Command: %s\n" % format_command(best.command))
    sys.stdout.write("Results: %s\n" % config.search.output_dir)
    return 0


def _create_backend(config: Any) -> Any:
    backend_name = config.backend.name.lower()
    if backend_name == "deephyper":
        return DeepHyperBackend(config)
    if backend_name == "gptune":
        return GPTuneBackend(config)
    if backend_name == "ytopt":
        return YtoptBackend(config)
    raise ValueError("unsupported tune backend: %s" % config.backend.name)
