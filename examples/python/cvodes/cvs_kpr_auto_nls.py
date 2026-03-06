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
# Kvaerno-Prothero-Robinson (KPR) ODE test problem
#
#   [u]' = [ a  b ] [ (-1 + u^2 - r(t)) / (2u) ] + [ r'(t) / (2u) ]
#   [v]    [ b  a ] [ (-2 + v^2 - s(t)) / (2v) ]   [ s'(t) / (2v) ]
#
# with analytic solution
#
#   u(t) = sqrt(1 + r(t))
#   v(t) = sqrt(2 + s(t))
#
# where
#
#   r(t) = 0.5 * cos(t)
#   s(t) = cos(20 * t)
#
# This example uses CVODE with the SUNNonlinearSolver_Auto module,
# which automatically switches between modified Newton and
# fixed-point iteration.
# -----------------------------------------------------------------

import sys
import argparse
import math
from pathlib import Path

from sundials4py.core import *
from sundials4py.cvodes import *


class KPRODE:
    def __init__(self, a=-2.0, b=0.5, c=0.5, d=-1.0):
        self.a = a
        self.b = b
        self.c = c
        self.d = d

    @staticmethod
    def r(t):
        return 0.5 * math.cos(t)

    @staticmethod
    def rdot(t):
        return -0.5 * math.sin(t)

    @staticmethod
    def s(t):
        return math.cos(20.0 * t)

    @staticmethod
    def sdot(t):
        return -20.0 * math.sin(20.0 * t)

    @classmethod
    def true_sol(cls, t):
        u = math.sqrt(1.0 + cls.r(t))
        v = math.sqrt(2.0 + cls.s(t))
        return u, v

    def set_init_cond(self, yvec):
        u0, v0 = self.true_sol(0.0)
        y = N_VGetArrayPointer(yvec)
        y[0] = u0
        y[1] = v0
        return 0

    def f(self, t, yvec, ydotvec, _user_data):
        y = N_VGetArrayPointer(yvec)
        ydot = N_VGetArrayPointer(ydotvec)

        u = float(y[0])
        v = float(y[1])

        tmp1 = (-1.0 + u * u - self.r(t)) / (2.0 * u)
        tmp2 = (-2.0 + v * v - self.s(t)) / (2.0 * v)

        ydot[0] = self.a * tmp1 + self.b * tmp2 + self.rdot(t) / (2.0 * u)
        ydot[1] = self.c * tmp1 + self.d * tmp2 + self.sdot(t) / (2.0 * v)
        return 0

    def jac(self, t, yvec, _fyvec, J, _user_data, _tmp1, _tmp2, _tmp3):
        y = N_VGetArrayPointer(yvec)
        u = float(y[0])
        v = float(y[1])

        Jdata = SUNDenseMatrix_Data(J)
        Jdata[0, 0] = self.a / 2.0 + (self.a * (1.0 + self.r(t)) - self.rdot(t)) / (2.0 * u * u)
        Jdata[1, 0] = self.c / 2.0 + (self.c * (1.0 + self.r(t))) / (2.0 * u * u)
        Jdata[0, 1] = self.b / 2.0 + (self.b * (2.0 + self.s(t))) / (2.0 * v * v)
        Jdata[1, 1] = self.d / 2.0 + (self.d * (2.0 + self.s(t)) - self.sdot(t)) / (2.0 * v * v)
        return 0


def parse_args(argv):

    parser = argparse.ArgumentParser(
        prog=argv[0],
        description="Kvaerno-Prothero-Robinson (KPR) ODE test problem using sundials4py.cvodes",
    )

    parser.add_argument("--rtol", type=float, default=1.0e-3, help="Relative tolerance.")
    parser.add_argument("--atol", type=float, default=1.0e-5, help="Absolute tolerance.")
    parser.add_argument("--dtout", type=float, default=1.0, help="Output interval.")
    parser.add_argument("--nout", type=int, default=10, help="Number of outputs.")
    parser.add_argument("--aa-depth", type=int, default=0, help="Anderson acceleration depth for Auto solver.")
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
        default="cv_kpr_info.log",
        help="SUNLogger INFO log output file (used for plotting solver regions).",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Disable plot generation (enabled by default if matplotlib is available).",
    )
    parser.add_argument(
        "--plot-file",
        type=str,
        default="cv_kpr_solution.png",
        help="Plot output filename.",
    )
    parser.add_argument(
        "--stepsize-plot-file",
        type=str,
        default="cv_kpr_stepsize.png",
        help="Step size plot output filename (requires SUNLogger INFO log).",
    )
    parser.add_argument(
        "--show-plot",
        action="store_true",
        help="Show plot in a window (also writes --plot-file).",
    )
    args, sundials_argv = parser.parse_known_args(argv[1:])

    return args, sundials_argv


def make_plots(args):
    if not args.no_plot:
        try:
            import matplotlib.pyplot as plt
            from matplotlib.patches import Patch
        except Exception as e:
            print(f"\nSkipping plot: matplotlib import failed ({e}).")
            return

        # Import repo-local suntools for parsing SUNLogger output
        repo_root = Path(__file__).resolve().parents[3]
        sys.path.insert(0, str(repo_root / "tools"))
        try:
            from suntools import logs as sunlog
        except Exception as e:
            print(f"\nSolver overlay disabled: failed to import suntools ({e}).")
            sunlog = None

        def solver_segments_from_log(log_path):
            if sunlog is None:
                return []
            if not Path(log_path).exists():
                return []
            try:
                log = sunlog.log_file_to_list(log_path)
            except Exception as e:
                print(f"\nSkipping solver overlay: failed to parse log file '{log_path}' ({e}).")
                return []

            segments = []
            for entry in log:
                if str(entry.get("status", "")) != "success":
                    continue
                tn = float(entry.get("tn", float("nan")))
                h = float(entry.get("h", float("nan")))
                nls = entry.get("nonlinear-solve", {})
                solver = nls.get("solver", None)
                if not solver or not math.isfinite(tn) or not math.isfinite(h):
                    continue
                segments.append((tn, tn + h, str(solver)))

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
            return merged

        segments = solver_segments_from_log(args.logfile)
        if not segments:
            for i in range(len(ts) - 1):
                d_fp = nfp_totals[i + 1] - nfp_totals[i]
                d_newt = nnewt_totals[i + 1] - nnewt_totals[i]
                if d_newt > d_fp:
                    solver = "Newton"
                elif d_fp > d_newt:
                    solver = "Fixed-Point"
                else:
                    solver = "Mixed"
                segments.append((ts[i], ts[i + 1], solver))

        colors = {"Newton": "tab:blue", "Fixed-Point": "tab:orange", "Mixed": "0.75"}

        fig, ax = plt.subplots(figsize=(9, 4.5))

        for t0, t1, solver in segments:
            c = colors.get(solver, "0.8")
            ax.axvspan(t0, t1, color=c, alpha=0.12, linewidth=0)

        ax.plot(ts, us, color="black", linewidth=1.8, label="u(t)")
        ax.plot(ts, vs, color="0.35", linewidth=1.8, label="v(t)")

        legend_patches = [
            Patch(facecolor=colors["Newton"], alpha=0.18, edgecolor="none", label="Newton steps"),
            Patch(facecolor=colors["Fixed-Point"], alpha=0.18, edgecolor="none", label="Fixed-Point steps"),
        ]
        if any(s == "Mixed" for _, _, s in segments):
            legend_patches.append(
                Patch(facecolor=colors["Mixed"], alpha=0.18, edgecolor="none", label="Mixed/unknown")
            )
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles + legend_patches, labels + [p.get_label() for p in legend_patches], loc="best")

        ax.set_xlabel("t")
        ax.set_ylabel("solution")
        ax.grid(alpha=0.3, linestyle="--")

        fig.tight_layout()
        plt.savefig(args.plot_file, bbox_inches="tight")

        # Step size vs step number plot (from SUNLogger info log)
        if sunlog is not None and Path(args.logfile).exists():
            log = sunlog.log_file_to_list(args.logfile)
            steps, _times, hs = sunlog.get_history(log, "h", "success")

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


def main(argv=None):
    if argv is None:
        argv = sys.argv

    args, sundials_argv = parse_args(argv)

    T0 = 0.0
    NEQ = 2
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

    coupling = 0.0
    stiffness = -1e2
    problem = KPRODE(a=stiffness, b=coupling, c=coupling, d=stiffness)
    problem.set_init_cond(y)

    cvode = CVodeCreate(CV_BDF, sunctx)
    assert cvode is not None

    status = CVodeInit(cvode.get(), problem.f, T0, y)
    assert status == CV_SUCCESS

    status = CVodeSStolerances(cvode.get(), args.rtol, args.atol)
    assert status == CV_SUCCESS

    status = CVodeSetMaxNumSteps(cvode.get(), 10000)
    assert status == CV_SUCCESS

    active_solver_type = (
        SUNNONLINSOL_AUTO_FIXEDPOINT if args.auto_init == "fixedpoint" else SUNNONLINSOL_AUTO_NEWTON
    )
    nls = SUNNonlinSol_Auto(y, args.aa_depth, active_solver_type, sunctx)
    assert nls is not None

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
    utrue0, vtrue0 = problem.true_sol(T0)

    print("\nKvaerno-Prothero-Robinson ODE test problem (sundials4py.cvodes):")
    print(f"    a = {problem.a}, b = {problem.b}, c = {problem.c}, d = {problem.d}")
    print(f"    nonlinear solver = SUNNonlinearSolver_Auto (init={args.auto_init}, aa_depth={args.aa_depth})")
    print(f"    reltol = {args.rtol:.2e}, abstol = {args.atol:.2e}\n")
    print("           t                   u                   v             |u - u*|            |v - v*|")
    print("   -------------------------------------------------------------------------------------------")
    print(
        f"  {T0:22.15e} {yarr[0]:22.15e} {yarr[1]:22.15e} "
        f"{abs(yarr[0] - utrue0):18.10e} {abs(yarr[1] - vtrue0):18.10e}"
    )

    ts = [T0]
    us = [float(yarr[0])]
    vs = [float(yarr[1])]
    _, nfp0, nnewt0 = SUNNonlinSolGetNumItersByType_Auto(nls)
    nfp_totals = [int(nfp0)]
    nnewt_totals = [int(nnewt0)]

    with open("cv_kpr_solution.txt", "w") as out:
        out.write("# t u v uerr verr\n")
        out.write(
            f"{T0:.16e} {yarr[0]:.16e} {yarr[1]:.16e} "
            f"{abs(yarr[0] - utrue0):.16e} {abs(yarr[1] - vtrue0):.16e}\n"
        )

        tout = T0 + args.dtout
        for _ in range(args.nout):
            status, tret = CVode(cvode.get(), tout, y, CV_NORMAL)
            assert status == CV_SUCCESS

            utrue, vtrue = problem.true_sol(tret)
            yarr = N_VGetArrayPointer(y)

            uerr = abs(yarr[0] - utrue)
            verr = abs(yarr[1] - vtrue)

            print(
                f"  {tret:22.15e} {yarr[0]:22.15e} {yarr[1]:22.15e} "
                f"{uerr:18.10e} {verr:18.10e}"
            )
            out.write(f"{tret:.16e} {yarr[0]:.16e} {yarr[1]:.16e} {uerr:.16e} {verr:.16e}\n")
            ts.append(float(tret))
            us.append(float(yarr[0]))
            vs.append(float(yarr[1]))
            _, nfp, nnewt = SUNNonlinSolGetNumItersByType_Auto(nls)
            nfp_totals.append(int(nfp))
            nnewt_totals.append(int(nnewt))

            tout = min(tout + args.dtout, Tf)

    print("   -------------------------------------------------------------------------------------------")

    status, nst = CVodeGetNumSteps(cvode.get())
    assert status == CV_SUCCESS
    status, nfe = CVodeGetNumRhsEvals(cvode.get())
    assert status == CV_SUCCESS
    status, nni = CVodeGetNumNonlinSolvIters(cvode.get())
    assert status == CV_SUCCESS
    status, ncfn = CVodeGetNumNonlinSolvConvFails(cvode.get())
    assert status == CV_SUCCESS
    status, nsetups = CVodeGetNumLinSolvSetups(cvode.get())
    assert status == CV_SUCCESS
    status, nje = CVodeGetNumJacEvals(cvode.get())
    assert status == CV_SUCCESS
    status, nfeLS = CVodeGetNumLinRhsEvals(cvode.get())
    assert status == CV_SUCCESS
    status, nfp, nnewt = SUNNonlinSolGetNumItersByType_Auto(nls)

    print("\nFinal Solver Statistics:")
    print(f"   Internal solver steps = {nst}")
    print(f"   Total RHS evals = {nfe}")
    print(f"   Total number of nonlinear solver iterations = {nni}")
    print(f"         newton = {nnewt}, fixed-point = {nfp}")
    print(f"   Total number of nonlinear solver convergence failures = {ncfn}")
    print(f"   Total number of Jacobian evaluations = {nje}")
    print(f"   Total linear solver setups = {nsetups}")
    print(f"   Total RHS evals for setting up the linear system = {nfeLS}")

    
def test_cvs_kpr_auto_nls():
    main(argv=["cvs_kpr_auto_nls.py", "--no-plot"])


if __name__ == "__main__":
    main()
