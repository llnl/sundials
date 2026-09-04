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

import asyncio
import csv
import json
import os
import re
import shlex
import statistics
import subprocess
import threading
import tempfile
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from suntools.tune.models import (
    ConstraintConfig,
    MetricConfig,
    ObjectiveConfig,
    ParameterSpec,
    TuneConfig,
)


@dataclass
class TrialResult:
    parameters: Dict[str, Any]
    metric: Optional[float]
    wall_time: float
    returncode: int
    stdout: str
    stderr: str
    command: List[str]
    success: bool
    error: Optional[str] = None
    constraint_metric: Optional[float] = None
    feasible: bool = True
    repetitions: int = 1


_ENV_VAR_PATTERN = re.compile(r"\$(\{[^}]+\}|[A-Za-z_][A-Za-z0-9_]*)")


def expand_environment_variables(value: str, environment: Dict[str, str]) -> str:
    """Expand shell-style environment variables without invoking a shell."""

    def replace(match: re.Match[str]) -> str:
        token = match.group(1)
        name = token[1:-1] if token.startswith("{") else token
        return environment.get(name, match.group(0))

    return os.path.expandvars(_ENV_VAR_PATTERN.sub(replace, value))


def build_trial_argv(
    command: str, args: Iterable[str], parameters: Iterable[ParameterSpec], values: Dict[str, Any]
) -> List[str]:
    argv = [command]
    argv.extend(str(arg) for arg in args)
    for parameter in parameters:
        argv.append(parameter.name)
        argv.extend(parameter.format_values(values[parameter.name]))
    return argv


def run_trial(config: TuneConfig, values: Dict[str, Any]) -> TrialResult:
    return _run_trial(config, values, config.parameters)


def run_baseline(config: TuneConfig) -> TrialResult:
    """Run the executable with its default settings and no tune parameters."""

    return _run_trial(config, {}, [])


def _run_trial(
    config: TuneConfig,
    values: Dict[str, Any],
    parameters: Iterable[ParameterSpec],
) -> TrialResult:
    env = os.environ.copy()
    env.update(config.executable.env)
    base_cwd = config.executable.cwd.resolve()
    command = _resolve_command_path(
        expand_environment_variables(config.executable.command, env), base_cwd
    )
    argv = build_trial_argv(
        command, config.executable.args, parameters, values
    )

    results = []
    with tempfile.TemporaryDirectory(prefix="suntools-trial-") as trial_dir:
        trial_cwd = Path(trial_dir)
        for _ in range(config.search.repetitions):
            start = time.perf_counter()
            proc = subprocess.run(
                argv,
                cwd=trial_dir,
                env=env,
                text=True,
                capture_output=True,
                check=False,
            )
            wall_time = time.perf_counter() - start
            results.append(
                _make_trial_result(
                    config,
                    values,
                    argv,
                    wall_time,
                    proc.returncode,
                    proc.stdout,
                    proc.stderr,
                    trial_cwd,
                )
            )

    return _average_trial_results(results)


def _resolve_command_path(command: str, base_cwd: Path) -> str:
    """Resolve path-like commands before running them from an isolated cwd."""

    command_path = Path(command)
    if command_path.is_absolute():
        return command
    if command_path.parent != Path(".") or command.startswith("."):
        return str((base_cwd / command_path).resolve())
    return command


def _average_trial_results(results: List[TrialResult]) -> TrialResult:
    """Average metrics over repeated executions of one sampled configuration."""

    successful = all(result.success for result in results)
    metric = (
        statistics.fmean(result.metric for result in results) if successful else None
    )
    constraint_metrics = [
        result.constraint_metric
        for result in results
        if result.constraint_metric is not None
    ]
    constraint_metric = (
        statistics.fmean(constraint_metrics)
        if len(constraint_metrics) == len(results)
        else None
    )
    errors = [
        "repetition %d: %s" % (index, result.error)
        for index, result in enumerate(results, start=1)
        if result.error
    ]
    return TrialResult(
        parameters=results[0].parameters,
        metric=metric,
        wall_time=statistics.fmean(result.wall_time for result in results),
        returncode=next(
            (result.returncode for result in results if result.returncode != 0), 0
        ),
        stdout="\n\n".join(result.stdout for result in results),
        stderr="\n\n".join(result.stderr for result in results),
        command=results[0].command,
        success=successful,
        error="; ".join(errors) if errors else None,
        constraint_metric=constraint_metric,
        feasible=successful and all(result.feasible for result in results),
        repetitions=len(results),
    )


async def run_trial_async(config: TuneConfig, values: Dict[str, Any]) -> TrialResult:
    """Run a trial without blocking an asyncio-based tuning evaluator."""

    result: List[TrialResult] = []
    error: List[BaseException] = []

    def run() -> None:
        try:
            result.append(run_trial(config, values))
        except BaseException as err:  # propagate trial startup failures
            error.append(err)

    thread = threading.Thread(target=run)
    thread.start()
    while thread.is_alive():
        await asyncio.sleep(0.01)
    thread.join()

    if error:
        raise error[0]
    return result[0]


def _make_trial_result(
    config: TuneConfig,
    values: Dict[str, Any],
    argv: List[str],
    wall_time: float,
    returncode: int,
    stdout: str,
    stderr: str,
    cwd: Optional[Path] = None,
) -> TrialResult:

    metric: Optional[float] = None
    constraint_metric: Optional[float] = None
    feasible = False
    error: Optional[str] = None
    metric_cwd = cwd or config.executable.cwd
    if returncode == 0:
        try:
            metric = extract_objective(
                config.objective, stdout, stderr, wall_time, metric_cwd
            )
        except ValueError as err:
            error = str(err)
        if metric is not None:
            feasible = config.constraint is None
        if metric is not None and config.constraint is not None:
            try:
                constraint_metric = extract_constraint(
                    config.constraint, stdout, stderr, metric_cwd
                )
                feasible = constraint_metric <= config.constraint.upper_bound
                if not feasible:
                    error = "constraint %s=%s exceeds upper bound %s" % (
                        config.constraint.metric,
                        constraint_metric,
                        config.constraint.upper_bound,
                    )
            except ValueError as err:
                error = str(err)
    else:
        error = "trial command exited with status %d" % returncode

    success = returncode == 0 and metric is not None
    return TrialResult(
        parameters=dict(values),
        metric=metric,
        wall_time=wall_time,
        returncode=returncode,
        stdout=stdout,
        stderr=stderr,
        command=argv,
        success=success,
        error=error,
        constraint_metric=constraint_metric,
        feasible=feasible,
    )


def extract_objective(
    objective: ObjectiveConfig,
    stdout: str,
    stderr: str,
    wall_time: float,
    cwd: Optional[Path] = None,
) -> float:
    return _extract_metric(objective, stdout, stderr, wall_time, cwd)


def extract_constraint(
    constraint: ConstraintConfig,
    stdout: str,
    stderr: str,
    cwd: Optional[Path] = None,
) -> float:
    return _extract_metric(constraint, stdout, stderr, None, cwd)


def _extract_metric(
    metric: MetricConfig,
    stdout: str,
    stderr: str,
    wall_time: Optional[float],
    cwd: Optional[Path],
) -> float:
    if not metric.regex:
        if metric.metric == "wall_time" and wall_time is not None:
            return wall_time
        raise ValueError("metric %s requires a regex" % metric.metric)

    source = metric.source or "stdout"
    if source == "stdout":
        text = stdout
    elif source == "stderr":
        text = stderr
    else:
        path = Path(source)
        if cwd is not None and not path.is_absolute():
            path = cwd / path
        with path.open("r") as fp:
            text = fp.read()

    patterns = [metric.regex] if isinstance(metric.regex, str) else metric.regex
    assert patterns is not None
    values = []
    for pattern in patterns:
        match = re.search(pattern, text, re.MULTILINE)
        if not match:
            raise ValueError("metric regex did not match %s" % source)
        try:
            raw_value = match.group(metric.group)
        except IndexError as err:
            raise ValueError("metric regex group was not found") from err
        values.append(float(raw_value))

    if metric.aggregation == "mean":
        return statistics.fmean(values)
    return sum(values)


def objective_to_score(
    direction: str, metric: Optional[float], feasible: bool = True
) -> float:
    if metric is None or not feasible:
        return float("-inf")
    if direction == "minimize":
        return -metric
    return metric


def select_best(results: Iterable[TrialResult], direction: str) -> Optional[TrialResult]:
    successful = [
        result
        for result in results
        if result.success and result.feasible and result.metric is not None
    ]
    if not successful:
        return None
    reverse = direction == "maximize"
    return sorted(successful, key=lambda result: result.metric, reverse=reverse)[0]


def select_worst(results: Iterable[TrialResult], direction: str) -> Optional[TrialResult]:
    successful = [
        result
        for result in results
        if result.success and result.feasible and result.metric is not None
    ]
    if not successful:
        return None
    reverse = direction == "minimize"
    return sorted(successful, key=lambda result: result.metric, reverse=reverse)[0]


def write_results(
    output_dir: Path,
    results: List[TrialResult],
    best: Optional[TrialResult],
    baseline: Optional[TrialResult] = None,
    worst: Optional[TrialResult] = None,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_trials_jsonl(output_dir / "trials.jsonl", results)
    _write_results_csv(output_dir / "results.csv", results)
    _write_best_json(output_dir / "best.json", best)
    _write_best_json(output_dir / "baseline.json", baseline)
    _write_best_json(output_dir / "worst.json", worst)


def _write_trials_jsonl(path: Path, results: List[TrialResult]) -> None:
    with path.open("w") as fp:
        for result in results:
            fp.write(json.dumps(asdict(result), sort_keys=True))
            fp.write("\n")


def _write_results_csv(path: Path, results: List[TrialResult]) -> None:
    parameter_names: List[str] = []
    for result in results:
        for name in result.parameters:
            if name not in parameter_names:
                parameter_names.append(name)

    fieldnames = [
        "trial",
        "success",
        "feasible",
        "metric",
        "constraint_metric",
        "repetitions",
        "wall_time",
        "returncode",
        "error",
    ] + parameter_names
    with path.open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=fieldnames)
        writer.writeheader()
        for index, result in enumerate(results):
            row: Dict[str, Any] = {
                "trial": index,
                "success": result.success,
                "feasible": result.feasible,
                "metric": "" if result.metric is None else result.metric,
                "constraint_metric": (
                    "" if result.constraint_metric is None else result.constraint_metric
                ),
                "repetitions": result.repetitions,
                "wall_time": result.wall_time,
                "returncode": result.returncode,
                "error": result.error or "",
            }
            row.update(result.parameters)
            writer.writerow(row)


def _write_best_json(path: Path, best: Optional[TrialResult]) -> None:
    payload: Dict[str, Any]
    if best is None:
        payload = {}
    else:
        payload = {
            "metric": best.metric,
            "constraint_metric": best.constraint_metric,
            "feasible": best.feasible,
            "repetitions": best.repetitions,
            "parameters": best.parameters,
            "command": best.command,
            "reproducible_command": format_command(best.command),
        }
    with path.open("w") as fp:
        json.dump(payload, fp, indent=2, sort_keys=True)
        fp.write("\n")


def format_command(argv: Iterable[str]) -> str:
    return " ".join(shlex.quote(str(arg)) for arg in argv)
