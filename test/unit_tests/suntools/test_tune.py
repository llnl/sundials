#!/usr/bin/env python3
# ---------------------------------------------------------------
# Programmer(s): Cody J. Balos
# ---------------------------------------------------------------
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
# ---------------------------------------------------------------

import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import textwrap
import types
import unittest
from unittest.mock import patch

from _testutils import add_repo_suntools_to_path


add_repo_suntools_to_path()

HAS_TUNE_DEPS = True
TUNE_SKIP_REASON = ""
try:
    from suntools.cli import build_parser
    from suntools.tune.config import config_from_args, load_config, parse_parameter_spec
    from suntools.tune.cli import _create_backend
    from suntools.tune.deephyper_backend import to_deephyper_problem
    from suntools.tune.gptune_backend import GPTuneBackend
    from suntools.tune.gptune_backend import to_gptune_problem
    from suntools.tune.models import (
        BackendConfig,
        ConstraintConfig,
        ExecutableConfig,
        ObjectiveConfig,
        ParameterSpec,
        TuneConfig,
    )
    from suntools.tune.runner import (
        TrialResult,
        build_trial_argv,
        expand_environment_variables,
        objective_to_score,
        run_baseline,
        run_trial,
        select_best,
        select_worst,
    )
    from suntools.tune.ytopt_backend import to_ytopt_problem
except ModuleNotFoundError as err:
    if err.name in ("pydantic", "yaml"):
        HAS_TUNE_DEPS = False
        TUNE_SKIP_REASON = "suntools tune dependencies are not installed"
    else:
        raise


def _write_script(directory, name, source):
    path = os.path.join(directory, name)
    with open(path, "w") as fp:
        fp.write(textwrap.dedent(source))
    return path


@unittest.skipUnless(HAS_TUNE_DEPS, TUNE_SKIP_REASON)
class TestTune(unittest.TestCase):
    def test_parse_float_log_parameter_spec(self):
        spec = parse_parameter_spec("cvode.nlscoef", "1e-4:0.3:log")
        self.assertEqual(spec.name, "cvode.nlscoef")
        self.assertEqual(spec.type, "float")
        self.assertEqual(spec.bounds, (1.0e-4, 0.3))
        self.assertEqual(spec.scale, "log")

    def test_parse_int_and_choice_parameter_specs(self):
        int_spec = parse_parameter_spec("arkode.order", "int:2:5")
        choice_spec = parse_parameter_spec("cvode.linsol", "choice:dense,spgmr")
        self.assertEqual(int_spec.type, "int")
        self.assertEqual(int_spec.bounds, (2.0, 5.0))
        self.assertEqual(choice_spec.type, "choice")
        self.assertEqual(choice_spec.values, ["dense", "spgmr"])

    def test_multi_value_choice_parameter_spec(self):
        spec = parse_parameter_spec(
            "arkode.table_names",
            "choice:ARKODE_DIRK_NONE ARKODE_ERK_NONE,"
            "ARKODE_DIRK_NONE ARKODE_ARK324L2SA_ERK_4_2_3,"
            "ARKODE_ARK324L2SA_DIRK_4_2_3 ARKODE_ARK324L2SA_ERK_4_2_3",
        )
        self.assertEqual(
            spec.values,
            [
                "ARKODE_DIRK_NONE ARKODE_ERK_NONE",
                "ARKODE_DIRK_NONE ARKODE_ARK324L2SA_ERK_4_2_3",
                "ARKODE_ARK324L2SA_DIRK_4_2_3 ARKODE_ARK324L2SA_ERK_4_2_3",
            ],
        )

    def test_invalid_parameter_spec_raises(self):
        with self.assertRaises(ValueError):
            parse_parameter_spec("cvode.nlscoef", "1e-4:0.3:badscale")
        with self.assertRaises(ValueError):
            ParameterSpec(name="cvode.nlscoef", type="float", bounds=(1.0, 0.5))

    def test_cli_args_build_tune_config(self):
        parser = build_parser()
        args = parser.parse_args(
            [
                "tune",
                "--params",
                "cvode.nlscoef",
                "1e-4:0.3:log",
                "--max-evals",
                "4",
                "--repetitions",
                "2",
                "--constraint-metric",
                "error",
                "--constraint-regex",
                r"error=([0-9.eE+-]+)",
                "--constraint-upper-bound",
                "1e-6",
                "./exe",
                "--rtol",
                "1e-6",
            ]
        )
        config = config_from_args(args)
        self.assertEqual(config.search.max_evals, 4)
        self.assertEqual(config.search.repetitions, 2)
        self.assertEqual(config.executable.command, "./exe")
        self.assertEqual(config.executable.args, ["--rtol", "1e-6"])
        self.assertEqual(config.parameters[0].name, "cvode.nlscoef")
        self.assertEqual(config.constraint.metric, "error")
        self.assertEqual(config.constraint.upper_bound, 1.0e-6)

    def test_yaml_validation(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = os.path.join(tmpdir, "tune.yaml")
            with open(config_path, "w") as fp:
                fp.write(
                    textwrap.dedent(
                        """
                        backend:
                          name: deephyper
                          options: {surrogate_model: ET, acq_func: UCBd}
                        search:
                          max_evals: 3
                          workers: 1
                          repetitions: 2
                          output_dir: tune-out
                        executable:
                          command: ./exe
                          args: ["--rtol", "1e-6"]
                          cwd: .
                          env: {}
                        parameters:
                          - name: cvode.nlscoef
                            type: float
                            bounds: [1.0e-4, 0.3]
                            scale: log
                        objective:
                          metric: wall_time
                          direction: minimize
                        constraint:
                          metric: error
                          source: stdout
                          regex: "error=([0-9.eE+-]+)"
                          upper_bound: 1.0e-6
                        """
                    )
                )
            config = load_config(config_path)
            self.assertEqual(config.backend.options["surrogate_model"], "ET")
            self.assertEqual(config.search.repetitions, 2)
            self.assertEqual(config.parameters[0].scale, "log")
            self.assertEqual(config.search.output_dir, Path(tmpdir) / "tune-out")
            self.assertIsInstance(config.constraint, ConstraintConfig)
            self.assertEqual(config.constraint.upper_bound, 1.0e-6)

    def test_build_trial_argv_appends_setoptions_pairs(self):
        parameters = [
            ParameterSpec(name="cvode.nlscoef", type="float", bounds=(1.0e-4, 0.3)),
            ParameterSpec(name="arkode.order", type="int", bounds=(2, 5)),
        ]
        argv = build_trial_argv(
            "./exe",
            ["--rtol", "1e-6"],
            parameters,
            {"cvode.nlscoef": 0.1, "arkode.order": 4.0},
        )
        self.assertEqual(
            argv,
            [
                "./exe",
                "--rtol",
                "1e-6",
                "cvode.nlscoef",
                "0.1",
                "arkode.order",
                "4",
            ],
        )

    def test_build_trial_argv_expands_multi_value_choice(self):
        parameters = [
            ParameterSpec(
                name="arkode.table_names",
                type="choice",
                values=[
                    "ARKODE_DIRK_NONE ARKODE_ERK_NONE",
                    "ARKODE_DIRK_NONE ARKODE_ARK324L2SA_ERK_4_2_3",
                    "ARKODE_ARK324L2SA_DIRK_4_2_3 ARKODE_ARK324L2SA_ERK_4_2_3",
                ],
            )
        ]
        argv = build_trial_argv(
            "./exe", [], parameters, {"arkode.table_names": "ARKODE_DIRK_NONE ARKODE_ERK_NONE"}
        )
        self.assertEqual(
            argv, ["./exe", "arkode.table_names", "ARKODE_DIRK_NONE", "ARKODE_ERK_NONE"]
        )

    def test_wall_time_objective(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            script = _write_script(
                tmpdir,
                "ok.py",
                """
                import sys
                sys.exit(0)
                """,
            )
            config = TuneConfig(
                executable=ExecutableConfig(
                    command=sys.executable, args=[script], cwd=tmpdir
                ),
                parameters=[
                    ParameterSpec(
                        name="cvode.nlscoef", type="float", bounds=(1.0e-4, 0.3)
                    )
                ],
            )
            result = run_trial(config, {"cvode.nlscoef": 0.1})
            self.assertTrue(result.success, result.error)
            self.assertIsNotNone(result.metric)
            self.assertGreaterEqual(result.metric, 0.0)

    def test_baseline_runs_without_tune_parameters(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            script = _write_script(
                tmpdir,
                "metric.py",
                """
                import sys
                print(f"objective={len(sys.argv[1:])}")
                """,
            )
            config = TuneConfig(
                search={"repetitions": 3},
                executable=ExecutableConfig(
                    command=sys.executable, args=[script], cwd=tmpdir
                ),
                parameters=[
                    ParameterSpec(name="x", type="float", bounds=(0.0, 1.0))
                ],
                objective=ObjectiveConfig(
                    metric="argument_count",
                    source="stdout",
                    regex=r"objective=([0-9]+)",
                ),
            )
            result = run_baseline(config)
            self.assertTrue(result.success, result.error)
            self.assertEqual(result.parameters, {})
            self.assertEqual(result.metric, 0.0)
            self.assertEqual(result.repetitions, 3)

    def test_repetitions_average_objective_and_constraint_metrics(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            script = _write_script(
                tmpdir,
                "metric.py",
                """
                from pathlib import Path
                counter = Path("counter")
                count = int(counter.read_text()) if counter.exists() else 0
                counter.write_text(str(count + 1))
                print(f"objective={count + 1}.0")
                print(f"error={0.1 + 0.02 * count}")
                """,
            )
            config = TuneConfig(
                search={"repetitions": 2},
                executable=ExecutableConfig(
                    command=sys.executable, args=[script], cwd=tmpdir
                ),
                parameters=[
                    ParameterSpec(name="x", type="float", bounds=(0.0, 1.0))
                ],
                objective=ObjectiveConfig(
                    metric="runtime",
                    source="stdout",
                    regex=r"objective=([0-9.eE+-]+)",
                ),
                constraint=ConstraintConfig(
                    metric="error",
                    source="stdout",
                    regex=r"error=([0-9.eE+-]+)",
                    upper_bound=0.2,
                ),
            )

            result = run_trial(config, {"x": 0.5})
            self.assertTrue(result.success, result.error)
            self.assertTrue(result.feasible)
            self.assertEqual(result.repetitions, 2)
            self.assertEqual(result.metric, 1.5)
            self.assertAlmostEqual(result.constraint_metric, 0.11)

    def test_command_path_expands_environment_variables(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            script = _write_script(
                tmpdir,
                "metric.py",
                """
                print("ok")
                """,
            )
            config = TuneConfig(
                executable=ExecutableConfig(
                    command="$PYTHON_BIN",
                    args=[script],
                    cwd=tmpdir,
                    env={"PYTHON_BIN": sys.executable, "TOOL_ROOT": tmpdir},
                ),
                parameters=[
                    ParameterSpec(
                        name="cvode.nlscoef", type="float", bounds=(1.0e-4, 0.3)
                    )
                ],
            )
            result = run_trial(config, {"cvode.nlscoef": 0.1})
            self.assertTrue(result.success, result.error)
            self.assertEqual(result.command[0], sys.executable)

            self.assertEqual(
                expand_environment_variables(
                    "${TOOL_ROOT}/metric.py", {"TOOL_ROOT": tmpdir}
                ),
                script,
            )

    def test_trial_uses_isolated_working_directory_for_file_metrics(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            script = _write_script(
                tmpdir,
                "metric.py",
                """
                from pathlib import Path
                Path("metric.txt").write_text("metric=2.5\\n")
                print(Path.cwd())
                """,
            )
            config = TuneConfig(
                executable=ExecutableConfig(
                    command=sys.executable, args=[script], cwd=tmpdir
                ),
                parameters=[
                    ParameterSpec(
                        name="cvode.nlscoef", type="float", bounds=(1.0e-4, 0.3)
                    )
                ],
                objective=ObjectiveConfig(
                    metric="file_metric",
                    source="metric.txt",
                    regex=r"metric=([0-9.eE+-]+)",
                ),
            )

            result = run_trial(config, {"cvode.nlscoef": 0.1})
            self.assertTrue(result.success, result.error)
            self.assertEqual(result.metric, 2.5)
            self.assertNotIn(str(Path(tmpdir)), result.stdout)
            self.assertFalse(Path(tmpdir, "metric.txt").exists())

    def test_parsed_stdout_metric(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            script = _write_script(
                tmpdir,
                "metric.py",
                """
                print("objective=3.25")
                """,
            )
            config = TuneConfig(
                executable=ExecutableConfig(
                    command=sys.executable, args=[script], cwd=tmpdir
                ),
                parameters=[
                    ParameterSpec(
                        name="cvode.nlscoef", type="float", bounds=(1.0e-4, 0.3)
                    )
                ],
                objective=ObjectiveConfig(
                    metric="residual",
                    direction="minimize",
                    source="stdout",
                    regex=r"objective=([0-9.]+)",
                    group=1,
                ),
            )
            result = run_trial(config, {"cvode.nlscoef": 0.1})
            self.assertTrue(result.success, result.error)
            self.assertEqual(result.metric, 3.25)

    def test_sum_of_multiple_stdout_metric_terms(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            script = _write_script(
                tmpdir,
                "metric.py",
                """
                print("Fe=2")
                print("Fi=5")
                print("nfeLS=3")
                print("nje=7")
                """,
            )
            config = TuneConfig(
                executable=ExecutableConfig(
                    command=sys.executable, args=[script], cwd=tmpdir
                ),
                parameters=[
                    ParameterSpec(
                        name="cvode.nlscoef", type="float", bounds=(1.0e-4, 0.3)
                    )
                ],
                objective=ObjectiveConfig(
                    metric="solver_work",
                    source="stdout",
                    regex=[
                        r"Fe=([0-9]+)",
                        r"Fi=([0-9]+)",
                        r"nfeLS=([0-9]+)",
                        r"nje=([0-9]+)",
                    ],
                    aggregation="sum",
                ),
            )
            result = run_trial(config, {"cvode.nlscoef": 0.1})
            self.assertTrue(result.success, result.error)
            self.assertEqual(result.metric, 17.0)

    def test_constraint_metric_filters_infeasible_trials(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            script = _write_script(
                tmpdir,
                "metric.py",
                """
                import sys
                values = dict(zip(sys.argv[1::2], sys.argv[2::2]))
                value = float(values["x"])
                print(f"objective={value}")
                print(f"error={2.0 * value}")
                """,
            )
            config = TuneConfig(
                executable=ExecutableConfig(
                    command=sys.executable, args=[script], cwd=tmpdir
                ),
                parameters=[
                    ParameterSpec(name="x", type="float", bounds=(0.0, 1.0))
                ],
                objective=ObjectiveConfig(
                    metric="runtime",
                    source="stdout",
                    regex=r"objective=([0-9.eE+-]+)",
                ),
                constraint=ConstraintConfig(
                    metric="error",
                    source="stdout",
                    regex=r"error=([0-9.eE+-]+)",
                    upper_bound=0.15,
                ),
            )

            infeasible = run_trial(config, {"x": 0.1})
            self.assertTrue(infeasible.success, infeasible.error)
            self.assertFalse(infeasible.feasible)
            self.assertEqual(infeasible.metric, 0.1)
            self.assertEqual(infeasible.constraint_metric, 0.2)

            feasible = run_trial(config, {"x": 0.05})
            self.assertTrue(feasible.success, feasible.error)
            self.assertTrue(feasible.feasible)
            self.assertEqual(feasible.constraint_metric, 0.1)

    def test_infeasible_results_are_excluded_from_best_and_score(self):
        infeasible = TrialResult(
            parameters={},
            metric=0.1,
            wall_time=0.1,
            returncode=0,
            stdout="",
            stderr="",
            command=[],
            success=True,
            feasible=False,
        )
        feasible = TrialResult(
            parameters={},
            metric=1.0,
            wall_time=1.0,
            returncode=0,
            stdout="",
            stderr="",
            command=[],
            success=True,
            feasible=True,
        )
        self.assertIs(select_best([infeasible, feasible], "minimize"), feasible)
        self.assertIs(select_worst([infeasible, feasible], "minimize"), feasible)
        self.assertEqual(objective_to_score("minimize", 0.1, False), float("-inf"))

    def test_failed_trial(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            script = _write_script(
                tmpdir,
                "fail.py",
                """
                import sys
                sys.exit(7)
                """,
            )
            config = TuneConfig(
                executable=ExecutableConfig(
                    command=sys.executable, args=[script], cwd=tmpdir
                ),
                parameters=[
                    ParameterSpec(
                        name="cvode.nlscoef", type="float", bounds=(1.0e-4, 0.3)
                    )
                ],
            )
            result = run_trial(config, {"cvode.nlscoef": 0.1})
            self.assertFalse(result.success)
            self.assertFalse(result.feasible)
            self.assertEqual(result.returncode, 7)
            self.assertIsNone(result.metric)

    def test_deephyper_problem_conversion_if_available(self):
        try:
            import deephyper  # noqa: F401
        except ModuleNotFoundError:
            self.skipTest("DeepHyper is not installed")

        problem = to_deephyper_problem(
            [
                ParameterSpec(
                    name="cvode.nlscoef",
                    type="float",
                    bounds=(1.0e-4, 0.3),
                    scale="log",
                ),
                ParameterSpec(name="arkode.order", type="int", bounds=(2, 5)),
                ParameterSpec(name="mode", type="choice", values=["a", "b"]),
            ]
        )
        self.assertIn("cvode.nlscoef", str(problem))
        self.assertIn("arkode.order", str(problem))
        self.assertIn("mode", str(problem))

    def test_ytopt_backend_can_be_selected(self):
        config = TuneConfig(
            backend=BackendConfig(name="ytopt"),
            executable=ExecutableConfig(command="./exe"),
            parameters=[
                ParameterSpec(name="cvode.nlscoef", type="float", bounds=(1.0e-4, 0.3))
            ],
        )
        self.assertEqual(_create_backend(config).__class__.__name__, "YtoptBackend")

    def test_gptune_backend_can_be_selected(self):
        config = TuneConfig(
            backend=BackendConfig(name="gptune"),
            executable=ExecutableConfig(command="./exe"),
            parameters=[
                ParameterSpec(name="cvode.nlscoef", type="float", bounds=(1.0e-4, 0.3))
            ],
        )
        self.assertEqual(_create_backend(config).__class__.__name__, "GPTuneBackend")

    def test_gptune_problem_conversion_with_problem_api(self):
        class FakeDimension:
            def __init__(self, *args, **kwargs):
                self.args = args
                self.kwargs = kwargs
                self.name = kwargs.get("name")

            def __repr__(self):
                return repr((self.args, self.kwargs))

        class FakeSpace:
            def __init__(self, dimensions):
                self.dimensions = dimensions

            def __iter__(self):
                return iter(self.dimensions)

            def __str__(self):
                return repr(self.dimensions)

        class FakeTuningProblem:
            def __init__(
                self,
                input_space,
                parameter_space,
                output_space,
                objective,
                constraints,
                models,
            ):
                self.input_space = input_space
                self.parameter_space = parameter_space
                self.output_space = output_space
                self.objective = objective
                self.constraints = constraints
                self.models = models

            def __str__(self):
                return repr(
                    (
                        self.input_space.dimensions,
                        self.parameter_space.dimensions,
                        self.output_space.dimensions,
                    )
                )

        fake_autotune = types.ModuleType("autotune")
        fake_problem = types.ModuleType("autotune.problem")
        fake_space = types.ModuleType("autotune.space")
        fake_gptune = types.ModuleType("GPTune")
        fake_data = types.ModuleType("GPTune.data")

        fake_problem.TuningProblem = FakeTuningProblem
        fake_space.Space = FakeSpace
        fake_space.Real = FakeDimension
        fake_space.Integer = FakeDimension
        fake_space.Categorical = FakeDimension
        fake_data.Categoricalnorm = FakeDimension

        modules = {
            "autotune": fake_autotune,
            "autotune.problem": fake_problem,
            "autotune.space": fake_space,
            "GPTune": fake_gptune,
            "GPTune.data": fake_data,
        }

        with patch.dict(sys.modules, modules):
            problem = to_gptune_problem(
                [
                    ParameterSpec(
                        name="cvode.nlscoef",
                        type="float",
                        bounds=(1.0e-4, 0.3),
                        scale="log",
                    ),
                    ParameterSpec(name="arkode.order", type="int", bounds=(2, 5)),
                    ParameterSpec(name="mode", type="choice", values=["a", "b"]),
                ]
            )
        self.assertIn("cvode.nlscoef", str(problem))
        self.assertIn("arkode.order", str(problem))
        self.assertIn("mode", str(problem))
        self.assertEqual(
            problem.parameter_space.dimensions[0].kwargs["prior"], "log-uniform"
        )
        self.assertEqual(problem.output_space.dimensions[0].kwargs["name"], "objective")

    def test_gptune_backend_run_with_fake_runtime(self):
        class FakeDimension:
            def __init__(self, *args, **kwargs):
                self.args = args
                self.kwargs = kwargs
                self.name = kwargs.get("name")

        class FakeSpace:
            def __init__(self, dimensions):
                self.dimensions = dimensions

        class FakeTuningProblem:
            def __init__(
                self,
                input_space,
                parameter_space,
                output_space,
                objective,
                constraints,
                models,
            ):
                self.input_space = input_space
                self.parameter_space = parameter_space
                self.output_space = output_space
                self.objective = objective
                self.constraints = constraints
                self.models = models

        class FakeComputer:
            def __init__(self, nodes=1, cores=1, hosts=None):
                self.nodes = nodes
                self.cores = cores
                self.hosts = hosts

        class FakeData:
            def __init__(self, problem):
                self.problem = problem

        class FakeOptions(dict):
            def validate(self, computer):
                self["validated"] = computer.cores

        class FakeGPTune:
            def __init__(self, problem, computer, data, options, **kwargs):
                self.problem = problem
                self.computer = computer
                self.data = data
                self.options = options
                self.kwargs = kwargs

            def MLA(self, NS, NS1, Tgiven, NI):
                samples = [0.2, 0.1]
                for index in range(NS):
                    self.problem.objective({"cvode.nlscoef": samples[index]})
                return self.data, None, {}

        fake_autotune = types.ModuleType("autotune")
        fake_problem = types.ModuleType("autotune.problem")
        fake_space = types.ModuleType("autotune.space")
        fake_gptune = types.ModuleType("GPTune")
        fake_computer = types.ModuleType("GPTune.computer")
        fake_data = types.ModuleType("GPTune.data")
        fake_gptune_mod = types.ModuleType("GPTune.gptune")
        fake_options = types.ModuleType("GPTune.options")

        fake_problem.TuningProblem = FakeTuningProblem
        fake_space.Space = FakeSpace
        fake_space.Real = FakeDimension
        fake_space.Integer = FakeDimension
        fake_space.Categorical = FakeDimension
        fake_computer.Computer = FakeComputer
        fake_data.Categoricalnorm = FakeDimension
        fake_data.Data = FakeData
        fake_gptune_mod.GPTune = FakeGPTune
        fake_options.Options = FakeOptions

        modules = {
            "autotune": fake_autotune,
            "autotune.problem": fake_problem,
            "autotune.space": fake_space,
            "GPTune": fake_gptune,
            "GPTune.computer": fake_computer,
            "GPTune.data": fake_data,
            "GPTune.gptune": fake_gptune_mod,
            "GPTune.options": fake_options,
        }

        with tempfile.TemporaryDirectory() as tmpdir, patch.dict(sys.modules, modules):
            script = _write_script(
                tmpdir,
                "metric.py",
                """
                import sys
                values = dict(zip(sys.argv[1::2], sys.argv[2::2]))
                print("objective=%s" % values.get("cvode.nlscoef", "0.3"))
                """,
            )
            output_dir = os.path.join(tmpdir, "out")
            config = TuneConfig(
                backend=BackendConfig(name="gptune"),
                search={"max_evals": 2, "workers": 3, "output_dir": output_dir},
                executable=ExecutableConfig(
                    command=sys.executable, args=[script], cwd=tmpdir
                ),
                parameters=[
                    ParameterSpec(
                        name="cvode.nlscoef", type="float", bounds=(1.0e-4, 0.3)
                    )
                ],
                objective=ObjectiveConfig(
                    metric="objective",
                    direction="minimize",
                    source="stdout",
                    regex=r"objective=([0-9.eE+-]+)",
                    group=1,
                ),
            )
            results = GPTuneBackend(config).run()

            self.assertEqual(len(results), 2)
            with open(os.path.join(output_dir, "best.json"), "r") as fp:
                best = json.load(fp)
            self.assertEqual(best["parameters"]["cvode.nlscoef"], 0.1)
            with open(os.path.join(output_dir, "baseline.json"), "r") as fp:
                baseline = json.load(fp)
            self.assertEqual(baseline["parameters"], {})
            self.assertEqual(baseline["metric"], 0.3)
            with open(os.path.join(output_dir, "worst.json"), "r") as fp:
                worst = json.load(fp)
            self.assertEqual(worst["parameters"]["cvode.nlscoef"], 0.2)

    def test_ytopt_problem_conversion_with_problem_api(self):
        class FakeProblem:
            def __init__(self):
                self.dimensions = []
                self.objectives = []

            def add_dim(self, name, domain):
                self.dimensions.append((name, domain))

            def add_objective(self, name):
                self.objectives.append(name)

            def __str__(self):
                return repr((self.dimensions, self.objectives))

        fake_ytopt = types.ModuleType("ytopt")
        fake_problem = types.ModuleType("ytopt.problem")
        fake_problem.Problem = FakeProblem

        modules = {
            "ytopt": fake_ytopt,
            "ytopt.problem": fake_problem,
        }

        with patch.dict(sys.modules, modules):
            problem = to_ytopt_problem(
                [
                    ParameterSpec(
                        name="cvode.nlscoef",
                        type="float",
                        bounds=(1.0e-4, 0.3),
                        scale="log",
                    ),
                    ParameterSpec(name="arkode.order", type="int", bounds=(2, 5)),
                    ParameterSpec(name="mode", type="choice", values=["a", "b"]),
                ]
            )
        self.assertIn("cvode.nlscoef", str(problem))
        self.assertIn("arkode.order", str(problem))
        self.assertIn("mode", str(problem))
        self.assertEqual(problem.objectives, ["objective"])

    def test_cli_integration_if_deephyper_available(self):
        try:
            import deephyper  # noqa: F401
        except ModuleNotFoundError:
            self.skipTest("DeepHyper is not installed")

        with tempfile.TemporaryDirectory() as tmpdir:
            script = _write_script(
                tmpdir,
                "fake_exe.py",
                """
                import sys
                values = dict(zip(sys.argv[1::2], sys.argv[2::2]))
                x = float(values.get("cvode.nlscoef", "0.3"))
                print(f"objective={x}")
                """,
            )
            output_dir = os.path.join(tmpdir, "out")
            env = os.environ.copy()
            env["PYTHONPATH"] = os.pathsep.join(sys.path)
            proc = subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "suntools.cli",
                    "tune",
                    "--params",
                    "cvode.nlscoef",
                    "1e-4:1e-3:log",
                    "--max-evals",
                    "2",
                    "--output-dir",
                    output_dir,
                    "--metric",
                    "objective",
                    "--objective-regex",
                    r"objective=([0-9.eE+-]+)",
                    sys.executable,
                    script,
                ],
                text=True,
                capture_output=True,
                check=False,
                env=env,
            )
            self.assertEqual(proc.returncode, 0, proc.stderr)
            with open(os.path.join(output_dir, "best.json"), "r") as fp:
                best = json.load(fp)
            self.assertIn("cvode.nlscoef", best["parameters"])


def run_tests():
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    suite.addTests(loader.loadTestsFromTestCase(TestTune))
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    return result.wasSuccessful()


if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)
