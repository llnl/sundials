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

import math
import os
import threading
from pathlib import Path
from typing import Any, Dict, List, Tuple

from suntools.tune.models import ParameterSpec, TuneConfig
from suntools.tune.runner import (
    TrialResult,
    run_baseline,
    run_trial,
    select_best,
    select_worst,
    write_results,
)


def to_gptune_problem(
    parameters: List[ParameterSpec], objective_name: str = "objective", objective=None
) -> Any:
    """Convert canonical parameter specs to a GPTune TuningProblem."""

    TuningProblem, Space, Real, Integer, Categorical = _import_gptune_problem_types()

    if objective is None:

        def objective(point):
            return [0.0]

    input_space = Space([Integer(0, 1, transform="normalize", name="suntools_task")])
    parameter_space = Space(
        [
            _make_gptune_parameter(parameter, Real, Integer, Categorical)
            for parameter in parameters
        ]
    )
    output_space = Space([Real(float("-inf"), float("inf"), name=objective_name)])
    return TuningProblem(input_space, parameter_space, output_space, objective, {}, None)


class GPTuneBackend:
    def __init__(self, config: TuneConfig):
        self.config = config
        self.baseline = None
        self.worst = None

    def run(self) -> List[TrialResult]:
        GPTune, Computer, Data, Options = _import_gptune_runtime()

        self.baseline = run_baseline(self.config)
        results: List[TrialResult] = []
        lock = threading.Lock()

        def objective(point: Dict[str, Any]) -> List[float]:
            values = _point_to_parameter_values(point, self.config.parameters)
            trial_result = run_trial(self.config, values)
            with lock:
                results.append(trial_result)
            return [
                _metric_to_gptune_objective(
                    self.config.objective.direction,
                    trial_result.metric,
                    trial_result.feasible,
                )
            ]

        problem = to_gptune_problem(
            self.config.parameters, self.config.objective.metric, objective
        )
        computer, options, constructor_options, run_options = _backend_options(
            Computer, Options, self.config.backend.options, self.config.search.workers
        )
        data = Data(problem)

        output_dir = Path(self.config.search.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        previous_cwd = os.getcwd()
        try:
            os.chdir(str(output_dir))
            tuner = GPTune(
                problem,
                computer=computer,
                data=data,
                options=options,
                **constructor_options,
            )
            _run_mla(tuner, self.config.search.max_evals, run_options)
        finally:
            os.chdir(previous_cwd)

        best = select_best(results, self.config.objective.direction)
        self.worst = select_worst(results, self.config.objective.direction)
        write_results(output_dir, results, best, self.baseline, self.worst)
        return results


def _import_gptune_problem_types() -> Tuple[Any, Any, Any, Any, Any]:
    try:
        from autotune.problem import TuningProblem
        from autotune.space import Categorical, Integer, Real, Space
    except ModuleNotFoundError as err:
        raise RuntimeError(
            "GPTune is required for the gptune tune backend. "
            "Install suntools with its project dependencies."
        ) from err

    try:
        from GPTune.data import Categoricalnorm

        Categorical = Categoricalnorm
    except ModuleNotFoundError:
        pass

    return TuningProblem, Space, Real, Integer, Categorical


def _import_gptune_runtime() -> Tuple[Any, Any, Any, Any]:
    try:
        from GPTune.computer import Computer
        from GPTune.data import Data
        from GPTune.gptune import GPTune
        from GPTune.options import Options
    except ModuleNotFoundError as err:
        raise RuntimeError(
            "GPTune is required for the gptune tune backend. "
            "Install suntools with its project dependencies."
        ) from err
    return GPTune, Computer, Data, Options


def _make_gptune_parameter(
    parameter: ParameterSpec, Real: Any, Integer: Any, Categorical: Any
) -> Any:
    if parameter.type == "choice":
        return Categorical(
            list(parameter.values or []), transform="onehot", name=parameter.name
        )

    assert parameter.bounds is not None
    lower, upper = parameter.bounds
    prior = "log-uniform" if parameter.scale == "log" else "uniform"
    if parameter.type == "int":
        return Integer(
            int(lower),
            int(upper),
            prior=prior,
            transform="normalize",
            name=parameter.name,
        )
    return Real(lower, upper, prior=prior, transform="normalize", name=parameter.name)


def _point_to_parameter_values(
    point: Dict[str, Any], parameters: List[ParameterSpec]
) -> Dict[str, Any]:
    values: Dict[str, Any] = {}
    for parameter in parameters:
        value = point[parameter.name]
        if parameter.type == "int":
            value = int(value)
        values[parameter.name] = value
    return values


def _metric_to_gptune_objective(
    direction: str, metric: Any, feasible: bool = True
) -> float:
    if metric is None or not feasible:
        return float("inf")
    value = float(metric)
    if direction == "maximize":
        return -value
    return value


def _backend_options(
    Computer: Any, Options: Any, backend_options: Dict[str, Any], workers: int
) -> Tuple[Any, Any, Dict[str, Any], Dict[str, Any]]:
    options = Options()
    option_values = dict(backend_options)
    nodes = int(option_values.pop("nodes", 1))
    cores = int(option_values.pop("cores", workers))
    hosts = option_values.pop("hosts", None)
    ns1 = option_values.pop("NS1", None)
    task = option_values.pop("task", [0])

    constructor_options: Dict[str, Any] = {
        "historydb": _as_bool(option_values.pop("historydb", False))
    }
    driverabspath = option_values.pop("driverabspath", None)
    if driverabspath is not None:
        constructor_options["driverabspath"] = driverabspath

    for name, value in option_values.items():
        options[name] = value
    _set_default_option(options, "distributed_memory_parallelism", False)
    _set_default_option(options, "shared_memory_parallelism", False)
    _set_default_option(options, "objective_evaluation_parallelism", workers > 1)
    _set_default_option(options, "objective_multisample_threads", workers)
    _set_default_option(options, "model_processes", 1)
    _set_default_option(options, "model_restarts", 1)
    _set_default_option(options, "verbose", False)

    computer = Computer(nodes=nodes, cores=cores, hosts=hosts)
    validate = getattr(options, "validate", None)
    if validate is not None:
        validate(computer=computer)

    run_options = {"Tgiven": [_task_to_list(task)], "NI": 1}
    if ns1 is not None:
        run_options["NS1"] = int(ns1)
    return computer, options, constructor_options, run_options


def _set_default_option(options: Any, name: str, value: Any) -> None:
    if name not in options or options[name] is None:
        options[name] = value


def _task_to_list(task: Any) -> List[Any]:
    if isinstance(task, (list, tuple)):
        return list(task)
    return [task]


def _as_bool(value: Any) -> Any:
    if isinstance(value, str):
        return value.lower() not in ("0", "false", "no", "off")
    return value


def _run_mla(tuner: Any, max_evals: int, run_options: Dict[str, Any]) -> None:
    ns1 = run_options.pop("NS1", max(1, min(max_evals, math.ceil(max_evals / 2))))
    tuner.MLA(NS=max_evals, NS1=ns1, **run_options)
