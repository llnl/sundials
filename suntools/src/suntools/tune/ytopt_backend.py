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
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Tuple

from suntools.tune.models import ParameterSpec, TuneConfig
from suntools.tune.runner import (
    TrialResult,
    objective_to_score,
    run_trial,
    select_best,
    write_results,
)


def to_ytopt_problem(
    parameters: List[ParameterSpec], objective_name: str = "objective"
) -> Any:
    """Convert canonical parameter specs to a Ytopt Problem."""

    try:
        from ytopt.problem import Problem
    except ModuleNotFoundError as err:
        raise RuntimeError(
            "Ytopt is required for the ytopt tune backend. "
            'Install suntools with the "ytopt" extra: python -m pip install "suntools[ytopt]".'
        ) from err

    problem = Problem()
    for parameter in parameters:
        _add_parameter(problem, parameter)
    _add_objective(problem, objective_name)
    return problem


class YtoptBackend:
    def __init__(self, config: TuneConfig):
        self.config = config

    def run(self) -> List[TrialResult]:
        AMBS, Evaluator = _import_ytopt()

        problem = to_ytopt_problem(self.config.parameters, self.config.objective.metric)
        results: List[TrialResult] = []
        lock = threading.Lock()

        def objective(sampled_values: Any) -> float:
            values = _sampled_values_to_dict(sampled_values)
            trial_result = run_trial(self.config, values)
            with lock:
                results.append(trial_result)
            return objective_to_score(
                self.config.objective.direction, trial_result.metric
            )

        _attach_objective(problem, objective)
        evaluator = _create_evaluator(
            Evaluator, problem, objective, self.config.search.workers
        )

        output_dir = Path(self.config.search.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        search_options = dict(self.config.backend.options)
        if "max_evals" not in search_options and _supports_keyword(AMBS, "max_evals"):
            search_options["max_evals"] = self.config.search.max_evals

        search = _create_search(AMBS, problem, evaluator, search_options)
        _run_search(search, self.config.search.max_evals)

        best = select_best(results, self.config.objective.direction)
        write_results(output_dir, results, best)
        return results


def _import_ytopt() -> Tuple[Any, Any]:
    try:
        from ytopt.search.ambs import AMBS
    except ModuleNotFoundError as err:
        raise RuntimeError(
            "Ytopt is required for the ytopt tune backend. "
            'Install suntools with the "ytopt" extra: python -m pip install "suntools[ytopt]".'
        ) from err

    try:
        from ytopt.evaluator.evaluate import Evaluator
    except ModuleNotFoundError:
        try:
            from ytopt.evaluator import Evaluator
        except ModuleNotFoundError as err:
            raise RuntimeError(
                "Ytopt is required for the ytopt tune backend. "
                'Install suntools with the "ytopt" extra: python -m pip install "suntools[ytopt]".'
            ) from err

    return AMBS, Evaluator


def _add_parameter(problem: Any, parameter: ParameterSpec) -> None:
    hyperparameter = _try_make_configspace_hyperparameter(parameter)
    if hyperparameter is not None:
        if _try_add_hyperparameter(problem, parameter, hyperparameter):
            return

    domain = _parameter_domain(parameter)
    if _try_add_domain(problem, parameter.name, domain):
        return

    input_space = getattr(problem, "input_space", None)
    if input_space is not None and _try_add_domain(input_space, parameter.name, domain):
        return

    raise RuntimeError("could not add ytopt parameter: %s" % parameter.name)


def _try_make_configspace_hyperparameter(parameter: ParameterSpec) -> Optional[Any]:
    try:
        return _make_configspace_hyperparameter(parameter)
    except ModuleNotFoundError:
        return None


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


def _try_add_hyperparameter(
    problem: Any, parameter: ParameterSpec, hyperparameter: Any
) -> bool:
    candidates = [
        problem,
        getattr(problem, "input_space", None),
        getattr(problem, "space", None),
    ]
    for target in candidates:
        if target is None:
            continue
        for method_name, args in (
            ("add_hyperparameter", (hyperparameter,)),
            ("add", (hyperparameter,)),
            ("add_dim", (parameter.name, hyperparameter)),
        ):
            method = getattr(target, method_name, None)
            if method is None:
                continue
            try:
                method(*args)
                return True
            except Exception:
                continue
    return False


def _parameter_domain(parameter: ParameterSpec) -> Any:
    if parameter.type == "choice":
        return list(parameter.values or [])

    assert parameter.bounds is not None
    lower, upper = parameter.bounds
    if parameter.type == "int":
        return (int(lower), int(upper))
    return (lower, upper)


def _try_add_domain(target: Any, name: str, domain: Any) -> bool:
    for method_name in ("add_dim", "add_parameter", "add_hyperparameter"):
        method = getattr(target, method_name, None)
        if method is None:
            continue
        for args in ((name, domain), (domain, name)):
            try:
                method(*args)
                return True
            except Exception:
                continue
    return False


def _add_objective(problem: Any, objective_name: str) -> None:
    candidates = [problem, getattr(problem, "output_space", None)]
    for target in candidates:
        if target is None:
            continue
        for method_name, args in (
            ("add_objective", (objective_name,)),
            ("add_dim", (objective_name, "float")),
            ("add", (objective_name,)),
        ):
            method = getattr(target, method_name, None)
            if method is None:
                continue
            try:
                method(*args)
                return
            except Exception:
                continue


def _attach_objective(problem: Any, objective: Callable[[Any], float]) -> None:
    for name in ("objective", "run", "evaluate"):
        try:
            setattr(problem, name, objective)
        except Exception:
            pass


def _sampled_values_to_dict(sampled_values: Any) -> Dict[str, Any]:
    if hasattr(sampled_values, "parameters"):
        sampled_values = sampled_values.parameters
    if hasattr(sampled_values, "get_dictionary"):
        sampled_values = sampled_values.get_dictionary()
    if isinstance(sampled_values, Mapping):
        return dict(sampled_values)
    return dict(sampled_values)


def _create_evaluator(
    Evaluator: Any,
    problem: Any,
    objective: Callable[[Any], float],
    workers: int,
) -> Any:
    create = getattr(Evaluator, "create", None)
    method_kwargs = {"num_workers": workers}
    attempts: List[Tuple[Tuple[Any, ...], Dict[str, Any]]] = []
    for method in ("thread", "threadPool", "serial"):
        attempts.extend(
            [
                ((problem,), {"method": method, "method_kwargs": method_kwargs}),
                (
                    (problem, objective),
                    {"method": method, "method_kwargs": method_kwargs},
                ),
                ((objective,), {"method": method, "method_kwargs": method_kwargs}),
            ]
        )
    attempts.extend([((problem,), {}), ((problem, objective), {}), ((objective,), {})])

    if create is not None:
        evaluator = _try_call(create, attempts)
        if evaluator is not None:
            return evaluator

    evaluator = _try_call(Evaluator, attempts)
    if evaluator is not None:
        return evaluator

    raise RuntimeError("could not create ytopt evaluator")


def _create_search(
    AMBS: Any, problem: Any, evaluator: Any, search_options: Dict[str, Any]
) -> Any:
    attempts = [
        ((), {"problem": problem, "evaluator": evaluator, **search_options}),
        ((problem, evaluator), search_options),
        ((), {"problem": problem, **search_options}),
        ((problem,), search_options),
    ]
    search = _try_call(AMBS, attempts)
    if search is None:
        raise RuntimeError("could not create ytopt AMBS search")
    return search


def _run_search(search: Any, max_evals: int) -> None:
    for method_name in ("search", "main", "run"):
        method = getattr(search, method_name, None)
        if method is None:
            continue
        attempts = [
            ((), {"max_evals": max_evals}),
            ((max_evals,), {}),
            ((), {}),
        ]
        if _try_call(method, attempts, allow_none=True) is not _CALL_FAILED:
            return
    raise RuntimeError("ytopt AMBS search object has no runnable search method")


_CALL_FAILED = object()


def _try_call(
    function: Callable[..., Any],
    attempts: Iterable[Tuple[Tuple[Any, ...], Dict[str, Any]]],
    allow_none: bool = False,
) -> Any:
    for args, kwargs in attempts:
        try:
            result = function(*args, **kwargs)
            if result is None and not allow_none:
                continue
            return result
        except TypeError:
            continue
    return _CALL_FAILED if allow_none else None


def _supports_keyword(function: Any, keyword: str) -> bool:
    try:
        signature = inspect.signature(function)
    except (TypeError, ValueError):
        return True
    return keyword in signature.parameters or any(
        parameter.kind == inspect.Parameter.VAR_KEYWORD
        for parameter in signature.parameters.values()
    )
