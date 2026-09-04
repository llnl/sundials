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

    baseline = getattr(backend, "baseline", None)
    if baseline is not None:
        _report_result("Baseline", baseline, config)

    best = select_best(results, config.objective.direction)
    if best is None:
        sys.stdout.write("Results: %s\n" % config.search.output_dir)
        if config.constraint is None:
            sys.stderr.write("error: no successful tune trials completed\n")
        else:
            sys.stderr.write("error: no feasible successful tune trials completed\n")
        return 1

    _report_result("Best", best, config)
    worst = getattr(backend, "worst", None)
    if worst is not None:
        _report_result("Worst", worst, config)
    sys.stdout.write("Results: %s\n" % config.search.output_dir)
    return 0


def _report_result(label: str, result: Any, config: Any) -> None:
    if result.metric is None:
        sys.stdout.write(
            "%s %s: failed (%s)\n"
            % (label, config.objective.metric, result.error or "metric unavailable")
        )
    else:
        sys.stdout.write("%s %s: %s\n" % (label, config.objective.metric, result.metric))
    if config.constraint is not None and result.constraint_metric is not None:
        sys.stdout.write(
            "  Constraint %s: %s <= %s\n"
            % (
                config.constraint.metric,
                result.constraint_metric,
                config.constraint.upper_bound,
            )
        )
    for name in sorted(result.parameters):
        sys.stdout.write("  %s = %s\n" % (name, result.parameters[name]))
    sys.stdout.write("  Command: %s\n" % format_command(result.command))


def _create_backend(config: Any) -> Any:
    backend_name = config.backend.name.lower()
    if backend_name == "deephyper":
        return DeepHyperBackend(config)
    if backend_name == "gptune":
        return GPTuneBackend(config)
    if backend_name == "ytopt":
        return YtoptBackend(config)
    raise ValueError("unsupported tune backend: %s" % config.backend.name)
