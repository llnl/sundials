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

from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import yaml

from suntools.tune.models import (
    BackendConfig,
    ExecutableConfig,
    ObjectiveConfig,
    ParameterSpec,
    SearchConfig,
    TuneConfig,
)


def parse_parameter_spec(name: str, spec: str) -> ParameterSpec:
    """Parse a CLI parameter specification into the canonical model."""

    if spec.startswith("choice:"):
        values = spec[len("choice:") :].split(",")
        return ParameterSpec(name=name, type="choice", values=values)

    if spec.startswith("int:"):
        parts = spec.split(":")
        if len(parts) != 3:
            raise ValueError("int parameter specs must use int:LOW:HIGH")
        lower = int(parts[1])
        upper = int(parts[2])
        return ParameterSpec(name=name, type="int", bounds=(lower, upper))

    parts = spec.split(":")
    if len(parts) not in (2, 3):
        raise ValueError("float parameter specs must use LOW:HIGH or LOW:HIGH:log")
    scale = "linear"
    if len(parts) == 3:
        if parts[2] != "log":
            raise ValueError("the only supported float scale suffix is :log")
        scale = "log"
    lower = float(parts[0])
    upper = float(parts[1])
    return ParameterSpec(name=name, type="float", bounds=(lower, upper), scale=scale)


def parse_key_value(items: Optional[Iterable[str]], option_name: str) -> Dict[str, str]:
    result: Dict[str, str] = {}
    if not items:
        return result
    for item in items:
        key, separator, value = item.partition("=")
        if not separator or not key:
            raise ValueError("%s entries must use KEY=VALUE" % option_name)
        result[key] = value
    return result


def parse_regex_group(value: Any) -> Any:
    if isinstance(value, str) and value.isdigit():
        return int(value)
    return value


def load_config(path: str) -> TuneConfig:
    config_path = Path(path)
    with config_path.open("r") as fp:
        data = yaml.safe_load(fp)
    if data is None:
        raise ValueError("empty tune configuration")
    config = TuneConfig.model_validate(data)
    return _resolve_relative_paths(config, config_path.parent)


def _resolve_relative_paths(config: TuneConfig, base_dir: Path) -> TuneConfig:
    output_dir = config.search.output_dir
    executable_cwd = config.executable.cwd
    updates: Dict[str, Any] = {}
    if not output_dir.is_absolute():
        updates["search"] = config.search.model_copy(
            update={"output_dir": base_dir / output_dir}
        )
    if not executable_cwd.is_absolute():
        updates["executable"] = config.executable.model_copy(
            update={"cwd": base_dir / executable_cwd}
        )
    if updates:
        config = config.model_copy(update=updates)
    return config


def config_from_args(args: Any) -> TuneConfig:
    if getattr(args, "config", None):
        return load_config(args.config)

    executable: Sequence[str] = getattr(args, "executable", None) or []
    if executable and executable[0] == "--":
        executable = executable[1:]
    if not executable:
        raise ValueError("suntools tune requires an executable command")

    parameter_items: Optional[List[Tuple[str, str]]] = getattr(args, "params", None)
    if not parameter_items:
        raise ValueError("suntools tune requires at least one --params KEY SPEC")
    parameters = [parse_parameter_spec(name, spec) for name, spec in parameter_items]

    objective = ObjectiveConfig(
        metric=args.metric,
        direction=args.direction,
        source=args.objective_source,
        regex=args.objective_regex,
        group=parse_regex_group(args.objective_group),
    )

    return TuneConfig(
        backend=BackendConfig(
            name=args.backend,
            options=parse_key_value(args.backend_option, "--backend-option"),
        ),
        search=SearchConfig(
            max_evals=args.max_evals,
            workers=args.workers,
            output_dir=Path(args.output_dir),
        ),
        executable=ExecutableConfig(
            command=executable[0],
            args=list(executable[1:]),
            cwd=Path(args.cwd),
            env=parse_key_value(args.env, "--env"),
        ),
        parameters=parameters,
        objective=objective,
    )
