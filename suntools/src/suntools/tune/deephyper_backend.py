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

import inspect
import threading
from pathlib import Path
from typing import Any, Callable, Dict, List

from suntools.tune.models import ParameterSpec, TuneConfig
from suntools.tune.runner import (
    TrialResult,
    objective_to_score,
    run_trial_async,
    select_best,
    write_results,
)


def to_deephyper_problem(parameters: List[ParameterSpec]) -> Any:
    """Convert canonical parameter specs to DeepHyper's HpProblem."""

    try:
        from deephyper.hpo import HpProblem
    except ModuleNotFoundError as err:
        raise RuntimeError(
            "DeepHyper is required for the deephyper tune backend. "
            "Install suntools with its project dependencies."
        ) from err

    problem = HpProblem()
    for parameter in parameters:
        _add_hyperparameter(problem, parameter)
    return problem


def _add_hyperparameter(problem: Any, parameter: ParameterSpec) -> None:
    try:
        hyperparameter = _make_configspace_hyperparameter(parameter)
        problem.add_hyperparameter(hyperparameter)
        return
    except Exception:
        pass

    if parameter.type == "choice":
        problem.add_hyperparameter(parameter.values, parameter.name)
    else:
        assert parameter.bounds is not None
        lower, upper = parameter.bounds
        if parameter.type == "int":
            value = (int(lower), int(upper))
        else:
            value = (lower, upper)
        problem.add_hyperparameter(value, parameter.name)


def _make_configspace_hyperparameter(parameter: ParameterSpec) -> Any:
    try:
        from ConfigSpace.hyperparameters import (
            CategoricalHyperparameter,
            UniformFloatHyperparameter,
            UniformIntegerHyperparameter,
        )
    except ModuleNotFoundError:
        from configspace.hyperparameters import (  # type: ignore[no-redef]
            CategoricalHyperparameter,
            UniformFloatHyperparameter,
            UniformIntegerHyperparameter,
        )

    if parameter.type == "choice":
        return CategoricalHyperparameter(parameter.name, choices=parameter.values)

    assert parameter.bounds is not None
    lower, upper = parameter.bounds
    log = parameter.scale == "log"
    if parameter.type == "int":
        return UniformIntegerHyperparameter(
            parameter.name, lower=int(lower), upper=int(upper), log=log
        )
    return UniformFloatHyperparameter(parameter.name, lower=lower, upper=upper, log=log)


class DeepHyperBackend:
    def __init__(self, config: TuneConfig):
        self.config = config

    def run(self) -> List[TrialResult]:
        try:
            from deephyper.evaluator import Evaluator
            from deephyper.hpo import CBO
        except ModuleNotFoundError as err:
            raise RuntimeError(
                "DeepHyper is required for the deephyper tune backend. "
                "Install suntools with its project dependencies."
            ) from err

        problem = to_deephyper_problem(self.config.parameters)
        results: List[TrialResult] = []
        lock = threading.Lock()

        async def objective(sampled_values: Dict[str, Any]) -> float:
            # DeepHyper 0.13 passes a RunningJob to evaluator callbacks while
            # older versions pass the sampled parameter dictionary directly.
            if not isinstance(sampled_values, dict) and hasattr(sampled_values, "parameters"):
                sampled_values = sampled_values.parameters
            trial_result = await run_trial_async(self.config, sampled_values)
            with lock:
                results.append(trial_result)
            return objective_to_score(self.config.objective.direction, trial_result.metric)

        evaluator = _create_evaluator(Evaluator, objective, self.config.search.workers)
        output_dir = Path(self.config.search.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        search_options = dict(self.config.backend.options)
        if "log_dir" not in search_options:
            search_options["log_dir"] = str(output_dir / "deephyper")
        if "evaluator" in inspect.signature(CBO.search).parameters:
            search = CBO(problem, **search_options)
            search.search(evaluator, max_evals=self.config.search.max_evals)
        else:
            search = CBO(problem, evaluator, **search_options)
            search.search(max_evals=self.config.search.max_evals)

        best = select_best(results, self.config.objective.direction)
        write_results(output_dir, results, best)
        return results


def _create_evaluator(
    Evaluator: Any, objective: Callable[[Dict[str, Any]], float], workers: int
) -> Any:
    return Evaluator.create(objective, method="serial", method_kwargs={"num_workers": workers})
