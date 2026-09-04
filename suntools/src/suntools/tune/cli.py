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

    _report_header(config)
    baseline = getattr(backend, "baseline", None)
    if baseline is not None:
        _report_result("Baseline", baseline, config)

    best = select_best(results, config.objective.direction)
    if best is None:
        if config.constraint is None:
            sys.stderr.write("error: no successful tune trials completed\n")
        else:
            sys.stderr.write("error: no feasible successful tune trials completed\n")
        return 1

    _report_result("Best", best, config)
    worst = getattr(backend, "worst", None)
    if worst is not None:
        _report_result("Worst", worst, config)
    return 0


def _report_header(config: Any) -> None:
    sys.stdout.write("suntools tune\n")
    sys.stdout.write("=============\n\n")
    sys.stdout.write(
        "Objective   : %s (%s)\n"
        % (config.objective.metric, config.objective.direction)
    )
    if config.constraint is not None:
        sys.stdout.write(
            "Constraint  : %s <= %s\n"
            % (config.constraint.metric, config.constraint.upper_bound)
        )
    sys.stdout.write(
        "Search      : %d evaluations, %d repetitions, %d worker%s\n"
        % (
            config.search.max_evals,
            config.search.repetitions,
            config.search.workers,
            "" if config.search.workers == 1 else "s",
        )
    )
    sys.stdout.write("Results     : %s\n\n" % config.search.output_dir)


def _report_result(label: str, result: Any, config: Any) -> None:
    sys.stdout.write("%s\n" % label)
    sys.stdout.write("%s\n" % ("-" * len(label)))
    sys.stdout.write("  Repetitions : %s\n" % result.repetitions)
    if result.metric is None:
        sys.stdout.write(
            "  Status      : failed (%s)\n"
            % (result.error or "metric unavailable")
        )
    else:
        sys.stdout.write("  Metric      : %s\n" % result.metric)
    if config.constraint is not None and result.constraint_metric is not None:
        sys.stdout.write(
            "  Constraint  : %s = %s (limit %s)\n"
            % (
                config.constraint.metric,
                result.constraint_metric,
                config.constraint.upper_bound,
            )
        )
    if result.parameters:
        sys.stdout.write("  Parameters  :\n")
        for name in sorted(result.parameters):
            sys.stdout.write("    %s = %s\n" % (name, result.parameters[name]))
    sys.stdout.write("  Command     : %s\n\n" % format_command(result.command))


def _create_backend(config: Any) -> Any:
    backend_name = config.backend.name.lower()
    if backend_name == "deephyper":
        return DeepHyperBackend(config)
    if backend_name == "gptune":
        return GPTuneBackend(config)
    if backend_name == "ytopt":
        return YtoptBackend(config)
    raise ValueError("unsupported tune backend: %s" % config.backend.name)
