#!/usr/bin/env python3
"""Compare warmed MPI heat solves, rotating backend order between samples."""

import argparse
import json
import math
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time

EXAMPLES = {
    "numpy": "ark_heat1D_mpi.py",
    "jax": "ark_heat1D_mpi_jax.py",
    "torch": "ark_heat1D_mpi_torch.py",
}


def validate(result, backend, ranks, nx, fields):
    if result["backend"] != backend or result["ranks"] != ranks:
        raise ValueError("Unexpected backend or rank count")
    if result["nx"] != nx or result["fields"] != fields:
        raise ValueError("Unexpected problem size")
    for key in ("solve_seconds", "max_error"):
        if not math.isfinite(result[key]) or result[key] < 0:
            raise ValueError(f"Invalid {key}")
    if result["error_limit"] != 2e-5:
        raise ValueError("Unexpected analytic error threshold")
    if result["status"] != 0 or result["max_error"] >= result["error_limit"]:
        raise ValueError("Solver failed or produced an incorrect solution")
    if not math.isclose(result["final_time"], result["tf"]):
        raise ValueError("Solver did not reach the final time")
    if result["dtype"] != "float64":
        raise ValueError("This comparison requires double precision")
    if backend != "numpy" and not result["device"].startswith("cuda"):
        raise ValueError("GPU backend did not use CUDA")
    if result["steps"] <= 0 or result["rhs_evals"] <= 0:
        raise ValueError("No solver work completed")
    if result["linear_failures"] or result["nonlinear_failures"]:
        raise ValueError("Solver convergence failures invalidate this comparison")


def verify_validation():
    """Keep negative controls so fast failures cannot masquerade as solves."""
    valid = dict(
        backend="jax",
        ranks=2,
        nx=201,
        fields=1,
        solve_seconds=1.0,
        max_error=1e-6,
        error_limit=2e-5,
        status=0,
        steps=10,
        rhs_evals=20,
        linear_failures=0,
        nonlinear_failures=0,
        final_time=0.1,
        tf=0.1,
        dtype="float64",
        device="cuda:0",
    )
    validate(valid, "jax", 2, 201, 1)
    for key, value in (
        ("backend", "numpy"),
        ("ranks", 1),
        ("nx", 99),
        ("fields", 2),
        ("solve_seconds", float("nan")),
        ("max_error", 1.0),
        ("error_limit", float("nan")),
        ("error_limit", float("inf")),
        ("max_error", 2e-5),
        ("status", -1),
        ("steps", 0),
        ("rhs_evals", 0),
        ("linear_failures", 1),
        ("nonlinear_failures", 1),
        ("final_time", 0),
        ("dtype", "float32"),
        ("device", "cpu"),
    ):
        try:
            validate(dict(valid, **{key: value}), "jax", 2, 201, 1)
        except ValueError:
            continue
        raise AssertionError(f"Validation failed to reject invalid {key}")


def main():
    verify_validation()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--ranks", type=int, nargs="+", default=[1, 2])
    parser.add_argument("--nx", type=int, default=201)
    parser.add_argument("--fields", type=int, nargs="+", default=[1, 256, 4096])
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--timeout", type=float, default=300)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if min(args.ranks + args.fields + [args.nx, args.repeats]) < 1:
        parser.error("Sizes and repeat count must be positive")
    directory = Path(__file__).resolve().parents[1] / "examples/python/arkode"
    env = os.environ.copy()
    env.update(
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
        XLA_PYTHON_CLIENT_PREALLOCATE="false",
        TRITON_DEFAULT_FP_FUSION="1",
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    backends = list(EXAMPLES)
    with args.output.open("w") as output:
        for ranks in args.ranks:
            for fields in args.fields:
                samples = {backend: [] for backend in backends}
                contract = None
                for repeat in range(args.repeats):
                    order = (
                        backends[repeat % len(backends) :]
                        + backends[: repeat % len(backends)]
                    )
                    for backend in order:
                        command = [
                            args.mpiexec,
                            "--bind-to",
                            "none",
                            "-n",
                            str(ranks),
                            args.python,
                            "-m",
                            "mpi4py",
                            str(directory / EXAMPLES[backend]),
                            "--nx",
                            str(args.nx),
                            "--fields",
                            str(fields),
                            "--json",
                            "--solver",
                            "sts",
                            "--multimode",
                        ]
                        start = time.perf_counter()
                        run = subprocess.run(
                            command,
                            env=env,
                            capture_output=True,
                            text=True,
                            timeout=args.timeout,
                        )
                        process_seconds = time.perf_counter() - start
                        if run.returncode:
                            raise RuntimeError(
                                f"{' '.join(command)}\n{run.stdout}\n{run.stderr}"
                            )
                        records = [
                            json.loads(line)
                            for line in run.stdout.splitlines()
                            if line.startswith("{")
                        ]
                        if len(records) != 1:
                            raise ValueError(
                                f"Expected one result, received: {run.stdout}"
                            )
                        result = records[0]
                        validate(result, backend, ranks, args.nx, fields)
                        specification = {
                            key: result[key]
                            for key in (
                                "nx",
                                "fields",
                                "rtol",
                                "atol",
                                "tf",
                                "k",
                                "dtype",
                                "solver",
                                "multimode",
                            )
                        }
                        if contract is not None and specification != contract:
                            raise ValueError(
                                "Backends used different numerical or solver settings"
                            )
                        contract = specification
                        result.update(
                            repeat=repeat,
                            process_seconds=process_seconds,
                            command=command,
                            stderr=run.stderr,
                        )
                        output.write(json.dumps(result) + "\n")
                        output.flush()
                        samples[backend].append(result["solve_seconds"])
                        print(
                            f"ranks={ranks} fields={fields} {backend}: {result['solve_seconds']:.6f} s",
                            flush=True,
                        )
                for backend, values in samples.items():
                    print(
                        f"ranks={ranks} fields={fields} {backend}: median={statistics.median(values):.6f} s "
                        f"range=[{min(values):.6f}, {max(values):.6f}] n={len(values)}",
                        flush=True,
                    )


if __name__ == "__main__":
    main()
