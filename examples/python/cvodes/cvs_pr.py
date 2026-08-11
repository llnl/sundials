#!/usr/bin/env python3
# -----------------------------------------------------------------
# Programmer(s): Cody J. Balos
# -----------------------------------------------------------------
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
# -----------------------------------------------------------------
# Scalar Prothero-Robinson (PR) ODE test problem
#
#   y' = a(t) (y - phi(t)) + phi'(t)
#
# with analytic solution
#
#   y(t) = phi(t)
#
# where
#
#   phi(t) = 1 + 0.5 * cos(t)
#
# This example uses CVODE with the SUNNonlinearSolver_Auto module,
# which automatically switches between modified Newton and
# fixed-point iteration.
# -----------------------------------------------------------------

import argparse
import math
import sys
import tempfile
from pathlib import Path

from sundials4py.core import *
from sundials4py.cvodes import *


def constant_stiffness(_t, base_stiffness, _t0, _tf, _args):
    return base_stiffness


def increasing_stiffness(t, base_stiffness, t0, tf, args):
    if tf <= t0:
        return base_stiffness

    tau = min(max((t - t0) / (tf - t0), 0.0), 1.0)
    scale = args.stiffness_growth_factor**tau
    return base_stiffness * scale


def oscillating_stiffness(t, base_stiffness, t0, tf, args):
    period = args.stiffness_period
    if period is None:
        period = max(tf - t0, sys.float_info.epsilon)

    phase = 2.0 * math.pi * (t - t0) / period
    scale = 10.0 ** (args.stiffness_oscillation_amplitude * math.sin(phase))
    return base_stiffness * scale


STIFFNESS_PROFILES = {
    "constant": constant_stiffness,
    "increasing": increasing_stiffness,
    "oscillating": oscillating_stiffness,
}


def make_stiffness_fn(args, t0, tf):
    stiffness_profile = STIFFNESS_PROFILES[args.stiffness_profile]

    def stiffness(t):
        return stiffness_profile(t, args.stiffness, t0, tf, args)

    return stiffness


class ProtheroRobinsonODE:
    def __init__(self, stiffness=-100.0, stiffness_fn=None):
        self.base_stiffness = stiffness
        self.stiffness_fn = stiffness_fn

    def stiffness(self, t):
        if self.stiffness_fn is None:
            return float(self.base_stiffness)
        return float(self.stiffness_fn(t))

    @staticmethod
    def phi(t):
        return 1.0 + 0.5 * math.cos(t)

    @staticmethod
    def phidot(t):
        return -0.5 * math.sin(t)

    @classmethod
    def true_sol(cls, t):
        return cls.phi(t)

    def set_init_cond(self, yvec):
        y = N_VGetArrayPointer(yvec)
        y[0] = self.true_sol(0.0)
        return 0

    def f(self, t, yvec, ydotvec, _user_data):
        y = N_VGetArrayPointer(yvec)
        ydot = N_VGetArrayPointer(ydotvec)

        ydot[0] = self.stiffness(t) * (float(y[0]) - self.phi(t)) + self.phidot(t)
        return 0

    def jac(self, t, _yvec, _fyvec, J, _user_data, _tmp1, _tmp2, _tmp3):
        Jdata = SUNDenseMatrix_Data(J)
        Jdata[0, 0] = self.stiffness(t)
        return 0


def parse_args(argv):

    parser = argparse.ArgumentParser(
        prog=argv[0],
        description="Scalar Prothero-Robinson (PR) ODE test problem using sundials4py.cvodes",
    )

    parser.add_argument("--stiffness", type=float, default=-1e2, help="Base stiffness parameter.")
    parser.add_argument(
        "--stiffness-profile",
        choices=sorted(STIFFNESS_PROFILES),
        default="constant",
        help=(
            "Time profile for a(t). 'constant' preserves the original PR problem; "
            "'increasing' grows |stiffness| geometrically over the integration interval; "
            "'oscillating' varies log10(|stiffness|) sinusoidally."
        ),
    )
    parser.add_argument(
        "--stiffness-growth-factor",
        type=float,
        default=100.0,
        help="Final/base stiffness-magnitude ratio for --stiffness-profile increasing.",
    )
    parser.add_argument(
        "--stiffness-oscillation-amplitude",
        type=float,
        default=0.75,
        help="Base-10 logarithmic amplitude for --stiffness-profile oscillating.",
    )
    parser.add_argument(
        "--stiffness-period",
        type=float,
        default=None,
        help="Oscillation period; defaults to the integration interval.",
    )
    parser.add_argument("--rtol", type=float, default=1.0e-3, help="Relative tolerance.")
    parser.add_argument("--atol", type=float, default=1.0e-5, help="Absolute tolerance.")
    parser.add_argument("--dtout", type=float, default=1.0, help="Output interval.")
    parser.add_argument("--nout", type=int, default=10, help="Number of outputs.")
    parser.add_argument(
        "--solver",
        choices=["auto", "newton", "fixedpoint"],
        default="auto",
        help="Nonlinear solver to use.",
    )
    parser.add_argument(
        "--aa-depth", type=int, default=0, help="Anderson acceleration depth for fixed-point paths."
    )
    parser.add_argument(
        "--auto-init",
        choices=["newton", "fixedpoint"],
        default="newton",
        help="Initial active solver type for SUNNonlinearSolver_Auto.",
    )
    parser.add_argument("--newt-to-fp-threshold", type=float, default=-1.0)
    parser.add_argument("--newt-to-fp-delay", type=int, default=-1)
    parser.add_argument("--fp-to-newt-threshold", type=float, default=-1.0)
    parser.add_argument("--fp-to-newt-delay", type=int, default=-1)
    parser.add_argument(
        "--logfile",
        type=str,
        default="cv_pr_info.log",
        help="SUNLogger INFO log output file (used for plotting solver regions).",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Disable plot generation (enabled by default if matplotlib is available).",
    )
    parser.add_argument(
        "--plot-file", type=str, default="cv_pr_solution.png", help="Plot output filename."
    )
    parser.add_argument(
        "--stepsize-plot-file",
        type=str,
        default="cv_pr_stepsize.png",
        help="Step size plot output filename (requires SUNLogger INFO log).",
    )
    parser.add_argument(
        "--stiffness-plot-file",
        type=str,
        default="cv_pr_stiffness.png",
        help="Stiffness plot output filename.",
    )
    parser.add_argument(
        "--show-plot", action="store_true", help="Show plot in a window (also writes --plot-file)."
    )
    args, sundials_argv = parser.parse_known_args(argv[1:])

    if args.stiffness_growth_factor <= 0.0:
        parser.error("--stiffness-growth-factor must be positive.")
    if args.stiffness_oscillation_amplitude < 0.0:
        parser.error("--stiffness-oscillation-amplitude must be nonnegative.")
    if args.stiffness_period is not None and args.stiffness_period <= 0.0:
        parser.error("--stiffness-period must be positive.")

    return args, sundials_argv


def _collect_solver_names(node, solvers, active_solvers):
    if isinstance(node, dict):
        solver = node.get("solver")
        active = node.get("active")
        if solver == "Auto" and active:
            active_solvers.add(str(active))
        elif isinstance(solver, str):
            solvers.add(solver)
        for value in node.values():
            _collect_solver_names(value, solvers, active_solvers)
    elif isinstance(node, list):
        for value in node:
            _collect_solver_names(value, solvers, active_solvers)


def _select_solver_name(solvers, active_solvers):
    solvers.discard("Auto")
    active_solvers.discard("Auto")

    for solver in ("Newton", "Fixed-Point"):
        if solver in active_solvers:
            return solver
    if len(solvers) == 1:
        return next(iter(solvers))
    if "Newton" in solvers:
        return "Newton"
    if "Fixed-Point" in solvers:
        return "Fixed-Point"
    return None


def _plot_solver_regions(ax, segments, colors):
    for t0, t1, solver in segments:
        c = colors.get(solver, "0.8")
        ax.axvspan(t0, t1, color=c, alpha=0.12, linewidth=0)


def _sample_stiffness(problem, t0, tf, nsamples=1001):
    if tf <= t0:
        return [t0], [abs(problem.stiffness(t0))]

    dt = (tf - t0) / (nsamples - 1)
    ts = [t0 + i * dt for i in range(nsamples)]
    return ts, [abs(problem.stiffness(t)) for t in ts]


def _load_solver_segments(args):
    try:
        from suntools import logs as sunlog
    except Exception as e:
        print(f"\nSolver overlay disabled: failed to import suntools ({e}).")
        return None, []

    if not Path(args.logfile).exists():
        return None, []

    try:
        log = sunlog.log_file_to_list(args.logfile)
    except Exception as e:
        print(f"\nSkipping solver overlay: failed to parse log file '{args.logfile}' ({e}).")
        return None, []

    segments = []
    for entry in log:
        if str(entry.get("status", "")) != "success":
            continue
        tn = float(entry.get("tn", float("nan")))
        h = float(entry.get("h", float("nan")))
        if not math.isfinite(tn) or not math.isfinite(h):
            continue

        solvers = set()
        active_solvers = set()
        _collect_solver_names(entry.get("nonlinear-solve", {}), solvers, active_solvers)
        solver = _select_solver_name(solvers, active_solvers)
        if solver is None:
            continue
        segments.append((tn, tn + h, solver))

    segments.sort(key=lambda x: x[0])
    merged = []
    for seg in segments:
        if not merged:
            merged.append(seg)
            continue
        a0, a1, aS = merged[-1]
        b0, b1, bS = seg
        if bS == aS and b0 <= a1:
            merged[-1] = (a0, max(a1, b1), aS)
        else:
            merged.append(seg)
    return log, merged


def make_plots(args, problem, ts, ys):
    if args.no_plot:
        return

    try:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch
    except Exception as e:
        print(f"\nSkipping plot: matplotlib import failed ({e}).")
        return

    log, segments = _load_solver_segments(args)
    colors = {"Newton": "tab:blue", "Fixed-Point": "tab:orange"}

    fig, ax = plt.subplots(figsize=(9, 4.5))
    _plot_solver_regions(ax, segments, colors)

    ytrue = [problem.true_sol(t) for t in ts]
    ax.plot(ts, ys, color="black", linewidth=1.8, label="y(t)")
    ax.plot(ts, ytrue, color="0.45", linestyle="--", linewidth=1.4, label="exact")

    legend_patches = [
        Patch(facecolor=colors["Newton"], alpha=0.18, edgecolor="none", label="Newton steps"),
        Patch(
            facecolor=colors["Fixed-Point"],
            alpha=0.18,
            edgecolor="none",
            label="Fixed-Point steps",
        ),
    ]
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles + legend_patches, labels + [p.get_label() for p in legend_patches])

    ax.set_xlabel("t")
    ax.set_ylabel("solution")
    ax.grid(alpha=0.3, linestyle="--")
    fig.tight_layout()
    plt.savefig(args.plot_file, bbox_inches="tight")

    fig3, ax3 = plt.subplots(figsize=(9, 4.0))
    _plot_solver_regions(ax3, segments, colors)
    stiffness_ts, stiffness_values = _sample_stiffness(problem, ts[0], ts[-1])
    ax3.plot(stiffness_ts, stiffness_values, color="tab:red", linewidth=1.8)
    ax3.set_xlabel("t")
    ax3.set_ylabel("|stiffness coefficient|")
    ax3.grid(alpha=0.3, linestyle="--")
    fig3.tight_layout()
    plt.savefig(args.stiffness_plot_file, bbox_inches="tight")

    if log is not None:
        try:
            from suntools import logs as sunlog

            steps, _times, hs = sunlog.get_history(log, "h", "success")
        except Exception as e:
            print(f"\nSkipping step size plot: failed to parse log file '{args.logfile}' ({e}).")
        else:
            fig2, ax2 = plt.subplots(figsize=(9, 4.0))
            ax2.plot(steps, hs, color="tab:green", marker=".", linewidth=1.2)
            ax2.set_xlabel("step")
            ax2.set_ylabel("step size (h)")
            ax2.grid(alpha=0.3, linestyle="--")
            fig2.tight_layout()
            plt.savefig(args.stepsize_plot_file, bbox_inches="tight")
    else:
        print("\nSkipping step size plot: requires SUNLogger INFO logfile and suntools.")

    if args.show_plot:
        plt.show()


def print_cvode_stats(cvode):
    with tempfile.NamedTemporaryFile(prefix="cv_pr_stats_", suffix=".txt", delete=False) as tmp:
        stats_path = Path(tmp.name)

    try:
        status, file_ptr = SUNFileOpen(str(stats_path), "w")
        assert status == CV_SUCCESS
        status = CVodePrintAllStats(cvode.get(), file_ptr, SUN_OUTPUTFORMAT_TABLE)
        assert status == CV_SUCCESS
        file_ptr = None
        print(stats_path.read_text(), end="")
    finally:
        stats_path.unlink(missing_ok=True)


def _create_nonlinear_solver(args, y, sunctx):
    if args.solver == "newton":
        return SUNNonlinSol_Newton(y, sunctx)

    if args.solver == "fixedpoint":
        return SUNNonlinSol_FixedPoint(y, args.aa_depth, sunctx)

    active_solver_type = (
        SUNNONLINSOL_AUTO_FIXEDPOINT
        if args.auto_init == "fixedpoint"
        else SUNNONLINSOL_AUTO_NEWTON
    )
    return SUNNonlinSol_Auto(y, args.aa_depth, active_solver_type, sunctx)


def main(argv=None):
    if argv is None:
        argv = sys.argv

    args, sundials_argv = parse_args(argv)

    T0 = 0.0
    NEQ = 1
    Tf = T0 + args.dtout * args.nout

    status, sunctx = SUNContext_Create(SUN_COMM_NULL)
    assert status == SUN_SUCCESS

    status, logger = SUNLogger_Create(SUN_COMM_NULL, 0)
    assert status == SUN_SUCCESS
    status = SUNLogger_SetInfoFilename(logger, args.logfile)
    assert status == SUN_SUCCESS
    status = SUNContext_SetLogger(sunctx, logger)
    assert status == SUN_SUCCESS

    y = N_VNew_Serial(NEQ, sunctx)
    assert y is not None

    stiffness_fn = make_stiffness_fn(args, T0, Tf)
    problem = ProtheroRobinsonODE(stiffness=args.stiffness, stiffness_fn=stiffness_fn)
    problem.set_init_cond(y)

    cvode = CVodeCreate(CV_BDF, sunctx)
    assert cvode is not None

    status = CVodeInit(cvode.get(), problem.f, T0, y)
    assert status == CV_SUCCESS

    status = CVodeSStolerances(cvode.get(), args.rtol, args.atol)
    assert status == CV_SUCCESS

    status = CVodeSetMaxNumSteps(cvode.get(), 10000)
    assert status == CV_SUCCESS

    nls = _create_nonlinear_solver(args, y, sunctx)
    assert nls is not None

    if args.solver == "auto":
        status = SUNNonlinSolSetSwitchingParameters_Auto(
            nls,
            args.newt_to_fp_threshold,
            args.newt_to_fp_delay,
            args.fp_to_newt_threshold,
            args.fp_to_newt_delay,
        )
        assert status == SUN_SUCCESS

    status = CVodeSetNonlinearSolver(cvode.get(), nls)
    assert status == CV_SUCCESS

    if args.solver in {"auto", "newton"}:
        A = SUNDenseMatrix(NEQ, NEQ, sunctx)
        LS = SUNLinSol_Dense(y, A, sunctx)
        assert A is not None and LS is not None

        status = CVodeSetLinearSolver(cvode.get(), LS, A)
        assert status == CV_SUCCESS

        status = CVodeSetJacFn(cvode.get(), problem.jac)
        assert status == CV_SUCCESS

    cvode_argv = [argv[0]] + sundials_argv
    status = CVodeSetOptions(cvode.get(), "", "", len(cvode_argv), cvode_argv)
    assert status == CV_SUCCESS

    yarr = N_VGetArrayPointer(y)
    ytrue0 = problem.true_sol(T0)

    print("\nProthero-Robinson ODE test problem (sundials4py.cvodes):")
    print(f"    stiffness profile = {args.stiffness_profile}")
    print(
        f"    a({T0:g}) = {problem.stiffness(T0):.6g}, "
        f"a({Tf:g}) = {problem.stiffness(Tf):.6g}"
    )
    print(f"    nonlinear solver = {args.solver}")
    if args.solver == "auto":
        print(f"    auto init = {args.auto_init}, aa_depth = {args.aa_depth}")
    print(f"    reltol = {args.rtol:.2e}, abstol = {args.atol:.2e}\n")
    print("           t                   y             |y - y*|")
    print("   ---------------------------------------------------------")
    print(f"  {T0:22.15e} {yarr[0]:22.15e} {abs(yarr[0] - ytrue0):18.10e}")

    ts = [T0]
    ys = [float(yarr[0])]
    with open("cv_pr_solution.txt", "w") as out:
        out.write("# t y yerr\n")
        out.write(f"{T0:.16e} {yarr[0]:.16e} {abs(yarr[0] - ytrue0):.16e}\n")

        tout = T0 + args.dtout
        for _ in range(args.nout):
            status, tret = CVode(cvode.get(), tout, y, CV_NORMAL)
            assert status == CV_SUCCESS

            ytrue = problem.true_sol(tret)
            yarr = N_VGetArrayPointer(y)
            yerr = abs(yarr[0] - ytrue)

            print(f"  {tret:22.15e} {yarr[0]:22.15e} {yerr:18.10e}")
            out.write(f"{tret:.16e} {yarr[0]:.16e} {yerr:.16e}\n")
            ts.append(float(tret))
            ys.append(float(yarr[0]))

            tout = min(tout + args.dtout, Tf)

    print("   ---------------------------------------------------------")

    print("\nFinal Solver Statistics:")
    print_cvode_stats(cvode)
    if args.solver == "auto":
        status, nfp, nnewt = SUNNonlinSolGetTotalNumItersByType_Auto(nls)
        assert status == CV_SUCCESS
        print(f"   Auto nonlinear solver iteration totals: newton = {nnewt}, fixed-point = {nfp}")

    make_plots(args, problem, ts, ys)


def test_cvs_pr_auto_nls():
    main(argv=["cvs_pr.py", "--no-plot"])


def test_cvs_pr_stiffness_profiles():
    args, _ = parse_args(["cvs_pr.py"])
    stiffness = make_stiffness_fn(args, 0.0, 10.0)
    assert math.isclose(stiffness(0.0), -100.0)
    assert math.isclose(stiffness(5.0), -100.0)
    assert math.isclose(stiffness(10.0), -100.0)

    args, _ = parse_args(
        [
            "cvs_pr.py",
            "--stiffness-profile",
            "increasing",
            "--stiffness",
            "-1.0",
            "--stiffness-growth-factor",
            "10.0",
        ]
    )
    stiffness = make_stiffness_fn(args, 0.0, 10.0)
    assert math.isclose(stiffness(0.0), -1.0)
    assert math.isclose(stiffness(5.0), -math.sqrt(10.0))
    assert math.isclose(stiffness(10.0), -10.0)

    args, _ = parse_args(
        [
            "cvs_pr.py",
            "--stiffness-profile",
            "oscillating",
            "--stiffness",
            "-10.0",
            "--stiffness-oscillation-amplitude",
            "0.5",
            "--stiffness-period",
            "4.0",
        ]
    )
    stiffness = make_stiffness_fn(args, 0.0, 10.0)
    assert math.isclose(stiffness(0.0), -10.0)
    assert math.isclose(stiffness(1.0), -10.0 * math.sqrt(10.0))
    assert math.isclose(stiffness(3.0), -10.0 / math.sqrt(10.0))


if __name__ == "__main__":
    main()
