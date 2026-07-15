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

import csv
import json
import os
import re
import shlex
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional

from suntools.tune.models import ObjectiveConfig, ParameterSpec, TuneConfig


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


def build_trial_argv(
    command: str,
    args: Iterable[str],
    parameters: Iterable[ParameterSpec],
    values: Dict[str, Any],
) -> List[str]:
    argv = [command]
    argv.extend(str(arg) for arg in args)
    for parameter in parameters:
        argv.append(parameter.name)
        argv.append(parameter.format_value(values[parameter.name]))
    return argv


def run_trial(config: TuneConfig, values: Dict[str, Any]) -> TrialResult:
    argv = build_trial_argv(
        config.executable.command,
        config.executable.args,
        config.parameters,
        values,
    )
    env = os.environ.copy()
    env.update(config.executable.env)

    start = time.perf_counter()
    proc = subprocess.run(
        argv,
        cwd=str(config.executable.cwd),
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    wall_time = time.perf_counter() - start

    metric: Optional[float] = None
    error: Optional[str] = None
    if proc.returncode == 0:
        try:
            metric = extract_objective(
                config.objective,
                proc.stdout,
                proc.stderr,
                wall_time,
                config.executable.cwd,
            )
        except ValueError as err:
            error = str(err)
    else:
        error = "trial command exited with status %d" % proc.returncode

    return TrialResult(
        parameters=dict(values),
        metric=metric,
        wall_time=wall_time,
        returncode=proc.returncode,
        stdout=proc.stdout,
        stderr=proc.stderr,
        command=argv,
        success=proc.returncode == 0 and metric is not None,
        error=error,
    )


def extract_objective(
    objective: ObjectiveConfig,
    stdout: str,
    stderr: str,
    wall_time: float,
    cwd: Optional[Path] = None,
) -> float:
    if not objective.regex:
        return wall_time

    source = objective.source or "stdout"
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

    match = re.search(objective.regex, text, re.MULTILINE)
    if not match:
        raise ValueError("objective regex did not match %s" % source)
    try:
        raw_value = match.group(objective.group)
    except IndexError as err:
        raise ValueError("objective regex group was not found") from err
    return float(raw_value)


def objective_to_score(direction: str, metric: Optional[float]) -> float:
    if metric is None:
        return float("-inf")
    if direction == "minimize":
        return -metric
    return metric


def select_best(results: Iterable[TrialResult], direction: str) -> Optional[TrialResult]:
    successful = [result for result in results if result.success and result.metric is not None]
    if not successful:
        return None
    reverse = direction == "maximize"
    return sorted(successful, key=lambda result: result.metric, reverse=reverse)[0]


def write_results(
    output_dir: Path, results: List[TrialResult], best: Optional[TrialResult]
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_trials_jsonl(output_dir / "trials.jsonl", results)
    _write_results_csv(output_dir / "results.csv", results)
    _write_best_json(output_dir / "best.json", best)


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
        "metric",
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
                "metric": "" if result.metric is None else result.metric,
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
            "parameters": best.parameters,
            "command": best.command,
            "reproducible_command": format_command(best.command),
        }
    with path.open("w") as fp:
        json.dump(payload, fp, indent=2, sort_keys=True)
        fp.write("\n")


def format_command(argv: Iterable[str]) -> str:
    return " ".join(shlex.quote(str(arg)) for arg in argv)
