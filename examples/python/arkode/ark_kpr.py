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
# This example uses ARKODE (ARKStep, implicit) with the
# SUNNonlinearSolver_Auto module, which automatically switches
# between modified Newton and fixed-point iteration.
# -----------------------------------------------------------------

import sys
import argparse
import math
from pathlib import Path

from sundials4py.core import *
from sundials4py.arkode import *


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

    def fi(self, t, yvec, ydotvec, _user_data):
        y = N_VGetArrayPointer(yvec)
        ydot = N_VGetArrayPointer(ydotvec)

        u = float(y[0])
        v = float(y[1])

        tmp1 = (-1.0 + u * u - self.r(t)) / (2.0 * u)
        tmp2 = (-2.0 + v * v - self.s(t)) / (2.0 * v)

        ydot[0] = self.a * tmp1 + self.b * tmp2 + self.rdot(t) / (2.0 * u)
        ydot[1] = self.c * tmp1 + self.d * tmp2 + self.sdot(t) / (2.0 * v)
        return 0

    def ji(self, t, yvec, _fyvec, J, _user_data, _tmp1, _tmp2, _tmp3):
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
        description="Kvaerno-Prothero-Robinson (KPR) ODE test problem using sundials4py.arkode (ARKStep) + SUNNonlinearSolver_Auto",
    )

    parser.add_argument("--rtol", type=float, default=1.0e-3, help="Relative tolerance.")
    parser.add_argument("--atol", type=float, default=1.0e-5, help="Absolute tolerance.")
    parser.add_argument("--dtout", type=float, default=1.0, help="Output interval.")
    parser.add_argument("--nout", type=int, default=10, help="Number of outputs.")
    parser.add_argument(
        "--aa-depth", type=int, default=0, help="Anderson acceleration depth for Auto solver."
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
        "--max-steps", type=int, default=10000, help="Maximum number of internal steps."
    )
    args, sundials_argv = parser.parse_known_args(argv[1:])
    return args, sundials_argv


def main(argv=None):
    if argv is None:
        argv = sys.argv

    args, sundials_argv = parse_args(argv)

    T0 = 0.0
    NEQ = 2
    Tf = T0 + args.dtout * args.nout

    status, sunctx = SUNContext_Create(SUN_COMM_NULL)
    assert status == SUN_SUCCESS

    y = N_VNew_Serial(NEQ, sunctx)
    assert y is not None

    coupling = 0.0
    stiffness = -1e2
    problem = KPRODE(a=stiffness, b=coupling, c=coupling, d=stiffness)
    problem.set_init_cond(y)

    # Create the ARKStep solver and configure tolerances/step limits.
    ark = ARKStepCreate(None, problem.fi, T0, y, sunctx)
    assert ark is not None

    status = ARKodeSStolerances(ark.get(), args.rtol, args.atol)
    assert status == ARK_SUCCESS

    status = ARKodeSetMaxNumSteps(ark.get(), args.max_steps)
    assert status == ARK_SUCCESS

    # Attach the switching nonlinear solver and the Newton linear solver pieces
    # needed whenever the auto method selects Newton.
    active_solver_type = (
        SUNNONLINSOL_AUTO_FIXEDPOINT
        if args.auto_init == "fixedpoint"
        else SUNNONLINSOL_AUTO_NEWTON
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

    status = ARKodeSetNonlinearSolver(ark.get(), nls)
    assert status == ARK_SUCCESS

    A = SUNDenseMatrix(NEQ, NEQ, sunctx)
    LS = SUNLinSol_Dense(y, A, sunctx)
    assert A is not None and LS is not None

    status = ARKodeSetLinearSolver(ark.get(), LS, A)
    assert status == ARK_SUCCESS

    status = ARKodeSetJacFn(ark.get(), problem.ji)
    assert status == ARK_SUCCESS

    arkode_argv = [argv[0]] + sundials_argv
    status = ARKodeSetOptions(ark.get(), "", "", len(arkode_argv), arkode_argv)
    assert status == ARK_SUCCESS

    yarr = N_VGetArrayPointer(y)
    utrue0, vtrue0 = problem.true_sol(T0)

    print("\nKvaerno-Prothero-Robinson ODE test problem (sundials4py.arkode, ARKStep implicit):")
    print(f"    a = {problem.a}, b = {problem.b}, c = {problem.c}, d = {problem.d}")
    print(
        f"    nonlinear solver = SUNNonlinearSolver_Auto (init={args.auto_init}, aa_depth={args.aa_depth})"
    )
    print(f"    reltol = {args.rtol:.2e}, abstol = {args.atol:.2e}\n")
    print(
        "           t                   u                   v             |u - u*|            |v - v*|"
    )
    print(
        "   -------------------------------------------------------------------------------------------"
    )
    print(
        f"  {T0:22.15e} {yarr[0]:22.15e} {yarr[1]:22.15e} "
        f"{abs(yarr[0] - utrue0):18.10e} {abs(yarr[1] - vtrue0):18.10e}"
    )

    ts = [T0]
    us = [float(yarr[0])]
    vs = [float(yarr[1])]
    with open("ark_kpr_solution.txt", "w") as out:
        out.write("# t u v uerr verr\n")
        out.write(
            f"{T0:.16e} {yarr[0]:.16e} {yarr[1]:.16e} "
            f"{abs(yarr[0] - utrue0):.16e} {abs(yarr[1] - vtrue0):.16e}\n"
        )

        tout = T0 + args.dtout
        for _ in range(args.nout):
            status, tret = ARKodeEvolve(ark.get(), tout, y, ARK_NORMAL)
            assert status == ARK_SUCCESS

            utrue, vtrue = problem.true_sol(tret)
            yarr = N_VGetArrayPointer(y)

            uerr = abs(yarr[0] - utrue)
            verr = abs(yarr[1] - vtrue)

            print(f"  {tret:22.15e} {yarr[0]:22.15e} {yarr[1]:22.15e} {uerr:18.10e} {verr:18.10e}")
            out.write(f"{tret:.16e} {yarr[0]:.16e} {yarr[1]:.16e} {uerr:.16e} {verr:.16e}\n")
            ts.append(float(tret))
            us.append(float(yarr[0]))
            vs.append(float(yarr[1]))

            tout = min(tout + args.dtout, Tf)

    print(
        "   -------------------------------------------------------------------------------------------"
    )

    print("\nFinal Solver Statistics:")
    status, file_ptr = SUNFileOpen("stdout", "w+")
    assert status == ARK_SUCCESS
    status = ARKodePrintAllStats(ark.get(), file_ptr, SUN_OUTPUTFORMAT_TABLE)
    assert status == ARK_SUCCESS

    status, nfp, nnewt = SUNNonlinSolGetTotalNumItersByType_Auto(nls)
    assert status == ARK_SUCCESS
    print(f"   Auto nonlinear solver iteration totals: newton = {nnewt}, fixed-point = {nfp}")


def test_ark_kpr_auto_nls():
    main()


if __name__ == "__main__":
    main()
