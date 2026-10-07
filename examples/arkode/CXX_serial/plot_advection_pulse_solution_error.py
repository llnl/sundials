#!/usr/bin/env python3
# -----------------------------------------------------------------------------
# Programmer(s): SUNDIALS development team
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
# Plot |u - u_ref| over time and space for one ark_advection_pulse run.
# ---------------------------------------------------------------------------

import argparse
import os
import subprocess
import tempfile
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "sundials-matplotlib")
)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


def default_executable():
    script_dir = Path(__file__).resolve().parent
    local_exe = script_dir / "ark_advection_pulse"
    if local_exe.exists():
        return local_exe

    repo_root = script_dir.parents[2]
    return repo_root / "build" / "examples" / "arkode" / "CXX_serial" / "ark_advection_pulse"


def run_solution(exe, workdir, args, solver, reltol, atol, ida_max_order):
    cmd = [
        str(exe),
        "--solver",
        solver,
        "--n",
        str(args.n),
        "--tf",
        f"{args.tf:.16g}",
        "--rtol",
        f"{reltol:.16g}",
        "--atol",
        f"{atol:.16g}",
        "--nout",
        str(args.nout),
        "--pulse-start",
        f"{args.pulse_start:.16g}",
        "--pulse-width",
        f"{args.pulse_width:.16g}",
    ]

    if solver == "pdae":
        cmd.extend(
            [
                "--partitions",
                str(args.partitions),
                "--pdae-h",
                f"{args.pdae_h:.16g}",
            ]
        )
    else:
        cmd.extend(["--ida-max-order", str(ida_max_order)])

    proc = subprocess.run(cmd, cwd=workdir, text=True, capture_output=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"{solver} run failed\n{proc.stdout}\n{proc.stderr}")

    path = Path(workdir) / "ark_advection_pulse.txt"
    if not path.exists():
        raise FileNotFoundError(f"{solver} run did not write {path}")

    return np.loadtxt(path)


def reshape_solution(data):
    times = np.unique(data[:, 0])
    xvals = np.unique(data[:, 1])
    expected = times.size * xvals.size
    if data.shape[0] != expected:
        raise ValueError("solution file does not contain a complete tensor grid")

    return times, xvals, data[:, 2].reshape((times.size, xvals.size))


def main():
    parser = argparse.ArgumentParser(
        description="Plot time-space error for one ark_advection_pulse run."
    )
    parser.add_argument("--exe", type=Path, default=default_executable())
    parser.add_argument("--solver", choices=["pdae", "ida"], default="pdae")
    parser.add_argument("--partitions", type=int, default=2)
    parser.add_argument("--pdae-h", type=float, default=1.0e-3)
    parser.add_argument("--n", type=int, default=201)
    parser.add_argument("--tf", type=float, default=1.2)
    parser.add_argument("--rtol", type=float, default=1.0e-6)
    parser.add_argument("--atol", type=float, default=1.0e-10)
    parser.add_argument("--ida-max-order", type=int, default=1)
    parser.add_argument("--ref-rtol", type=float, default=1.0e-10)
    parser.add_argument("--ref-atol", type=float, default=1.0e-12)
    parser.add_argument("--ref-ida-max-order", type=int, default=0)
    parser.add_argument("--nout", type=int, default=80)
    parser.add_argument("--pulse-start", type=float, default=0.05)
    parser.add_argument("--pulse-width", type=float, default=0.20)
    parser.add_argument("--error-floor", type=float, default=1.0e-16)
    parser.add_argument("--csv", type=Path, default=Path("ark_advection_pulse_solution_error.csv"))
    parser.add_argument("--png", type=Path, default=Path("ark_advection_pulse_solution_error.png"))
    args = parser.parse_args()

    exe = args.exe.resolve()
    if not exe.exists():
        raise FileNotFoundError(f"Could not find executable: {exe}")

    with tempfile.TemporaryDirectory(prefix="ark_advection_pulse_run_") as run_dir:
        run_data = run_solution(
            exe, run_dir, args, args.solver, args.rtol, args.atol, args.ida_max_order
        )

    with tempfile.TemporaryDirectory(prefix="ark_advection_pulse_ref_") as ref_dir:
        ref_data = run_solution(
            exe,
            ref_dir,
            args,
            "ida",
            args.ref_rtol,
            args.ref_atol,
            args.ref_ida_max_order,
        )

    times, xvals, u = reshape_solution(run_data)
    ref_times, ref_xvals, uref = reshape_solution(ref_data)

    if not np.allclose(times, ref_times) or not np.allclose(xvals, ref_xvals):
        raise ValueError("run and reference solution grids do not match")

    err = np.abs(u - uref)

    np.savetxt(
        args.csv,
        np.column_stack(
            (
                np.repeat(times, xvals.size),
                np.tile(xvals, times.size),
                err.reshape(times.size * xvals.size),
            )
        ),
        header="t x abs_error",
        comments="",
    )

    tgrid, xgrid = np.meshgrid(times, xvals, indexing="ij")
    zvals = np.log10(np.maximum(err, args.error_floor))

    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(tgrid, xgrid, zvals, cmap="viridis", linewidth=0, antialiased=True)
    ax.set_xlabel("t")
    ax.set_ylabel("x")
    ax.set_zlabel("log10(|u - u_ref|)")
    ax.set_title("Advection pulse solution error")
    ax.view_init(elev=28, azim=-135)
    fig.colorbar(surf, shrink=0.65, pad=0.12, label="log10(error)")
    fig.tight_layout()
    fig.savefig(args.png, dpi=200)

    print(f"max error = {np.max(err):.16e}")
    print(f"wrote {args.csv}")
    print(f"wrote {args.png}")


if __name__ == "__main__":
    main()
