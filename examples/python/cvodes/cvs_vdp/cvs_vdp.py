#!/usr/bin/env python3
# -----------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
# -----------------------------------------------------------------
# We solve the classic Van der Pol problem:
#   y'' - mu*(1 - y^2)*y' + y = 0,  y(0) = 2,  y'(0) = 0.
# This second-order ODE is converted to a first-order system by defining
#   y0 = y,  y1 = y'
# giving
#   y0' = y1
#   y1' = mu*(1 - y0^2)*y1 - y0.
# We use the SUNNonlinearSolver_Auto module to solve the implicit
# system. This solver automatically switches between modified Newton
# iteration and fixed-point iteration using a stiffness metric.
# -----------------------------------------------------------------

import sys
import argparse
import numpy as np
from sundials4py.core import *
from sundials4py.cvodes import *


class VanDerPolODE:
    def __init__(self, mu, y10, y20):
        self.mu = mu
        self.y10 = y10
        self.y20 = y20

    def set_init_cond(self, yvec):
        y = N_VGetArrayPointer(yvec)
        y[0] = self.y10
        y[1] = self.y20
        return 0

    def f(self, t, yvec, ydotvec):
        y = N_VGetArrayPointer(yvec)
        ydot = N_VGetArrayPointer(ydotvec)
        mu = self.mu
        ydot[0] = y[1]
        ydot[1] = mu * (1.0 - y[0] ** 2) * y[1] - y[0]
        return 0

    def jac(self, t, yvec, fyvec, J, tmp1, tmp2, tmp3):
        y = N_VGetArrayPointer(yvec)
        mu = self.mu
        Jdata = SUNDenseMatrix_Data(J)
        Jdata[0, 0] = 0.0
        Jdata[0, 1] = 1.0
        Jdata[1, 0] = -2.0 * mu * y[0] * y[1] - 1.0
        Jdata[1, 1] = mu * (1.0 - y[0] ** 2)
        return 0


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(
        prog=argv[0], description="Van der Pol oscillator example using sundials4py.cvodes."
    )
    parser.add_argument(
        "--mu", type=float, default=100.0, help="Van der Pol stiffness parameter (default: 100.0)"
    )
    parser.add_argument(
        "--Tf", "--tf", dest="tf", type=float, default=10.0, help="Final time (default: 10.0)"
    )
    parser.add_argument(
        "--solver",
        choices=["newton", "fixedpoint", "auto"],
        default="auto",
        help="Which nonlinear solver to use",
    )
    parser.add_argument(
        "--plot", action="store_true", help="Generate a plot of the solution (requires matplotlib)"
    )
    parser.add_argument(
        "--plot-file",
        default="cv_vdp_solution.png",
        help="Output filename for plot (default: cv_vdp_solution.png)",
    )
    parser.add_argument(
        "--show", action="store_true", help="Display plot window (in addition to saving)"
    )
    parser.set_defaults(newton=True)
    args, sundials_argv = parser.parse_known_args(argv[1:])

    mu = args.mu
    y10 = 2.0
    y20 = 0.0
    T0 = 0.0
    Tf = args.tf
    dTout = 1.0
    NEQ = 2
    Nt = int(np.ceil(Tf / dTout))
    reltol = 1e-4
    abstol = 1e-4
    aa_depth = 0

    status, sunctx = SUNContext_Create(SUN_COMM_NULL)
    y = N_VNew_Serial(NEQ, sunctx)

    ode = VanDerPolODE(mu, y10, y20)
    ode.set_init_cond(y)

    cvode = CVodeCreate(CV_BDF, sunctx)

    status = CVodeInit(cvode.get(), lambda t, y, ydot, _: ode.f(t, y, ydot), T0, y)
    assert status == CV_SUCCESS

    status = CVodeSStolerances(cvode.get(), reltol, abstol)
    assert status == CV_SUCCESS

    status = CVodeSetMaxNumSteps(cvode.get(), 10000)
    assert status == CV_SUCCESS

    if args.solver == "newton":
        NLS = SUNNonlinSol_Newton(y, sunctx)
    elif args.solver == "fixedpoint":
        NLS = SUNNonlinSol_FixedPoint(y, 0, sunctx)
    else:
        NLS = SUNNonlinSol_Auto(y, aa_depth, SUNNONLINSOL_AUTO_NEWTON, sunctx)

    status = CVodeSetNonlinearSolver(cvode.get(), NLS)
    assert status == CV_SUCCESS

    if args.solver == "newton" or args.solver == "auto":
        A = SUNDenseMatrix(NEQ, NEQ, sunctx)
        LS = SUNLinSol_Dense(y, A, sunctx)

        status = CVodeSetLinearSolver(cvode.get(), LS, A)
        assert status == CV_SUCCESS

        status = CVodeSetJacFn(
            cvode.get(),
            lambda t, yvec, fyvec, J, tmp1, tmp2, tmp3, _: ode.jac(
                t, yvec, fyvec, J, tmp1, tmp2, tmp3
            ),
        )
        assert status == CV_SUCCESS

    cvode_argv = [argv[0]] + sundials_argv
    status = CVodeSetOptions(cvode.get(), "", "", len(cvode_argv), cvode_argv)
    assert status == CV_SUCCESS

    yarr = N_VGetArrayPointer(y)
    print("\nVan der Pol oscillator solved with sundials4py.cvodes:")
    print(f"    initial conditions: y1 = {y10}, y2 = {y20}")
    print(f"    mu = {mu}")
    print(f"    nonlinear solver = {args.solver}")
    print(f"    reltol = {reltol}, abstol = {abstol}\n")
    print("        t           y1           y2")
    print("   -----------------------------------")
    print(f"  {T0:10.6f}  {yarr[0]:10.6f}  {yarr[1]:10.6f}")

    ts = [T0]
    y1s = [float(yarr[0])]
    y2s = [float(yarr[1])]

    with open("cv_vdp_solution.txt", "w") as UFID:
        UFID.write("# t y1 y2\n")
        UFID.write(f" {T0:.16e} {yarr[0]:.16e} {yarr[1]:.16e}\n")
        tout = T0 + dTout
        for iout in range(Nt):
            status, tret = CVode(cvode.get(), tout, y, CV_NORMAL)
            yarr = N_VGetArrayPointer(y)
            print(f"  {tret:10.6f}  {yarr[0]:10.6f}  {yarr[1]:10.6f}")
            UFID.write(f" {tret:.16e} {yarr[0]:.16e} {yarr[1]:.16e}\n")
            ts.append(float(tret))
            y1s.append(float(yarr[0]))
            y2s.append(float(yarr[1]))
            if status == CV_SUCCESS:
                tout += dTout
                tout = min(tout, Tf)
            else:
                print("Solver failure, stopping integration")
                break
        print("   -----------------------------------")

    status, nst = CVodeGetNumSteps(cvode.get())
    assert status == CV_SUCCESS
    status, nfe = CVodeGetNumRhsEvals(cvode.get())
    assert status == CV_SUCCESS
    status, nni = CVodeGetNumNonlinSolvIters(cvode.get())
    assert status == CV_SUCCESS
    status, ncfn = CVodeGetNumNonlinSolvConvFails(cvode.get())
    assert status == CV_SUCCESS

    if args.solver == "newton" or args.solver == "auto":
        status, nsetups = CVodeGetNumLinSolvSetups(cvode.get())
        assert status == CV_SUCCESS
        status, nje = CVodeGetNumJacEvals(cvode.get())
        assert status == CV_SUCCESS
        status, nfeLS = CVodeGetNumLinRhsEvals(cvode.get())
        assert status == CV_SUCCESS

    print("\nFinal Solver Statistics:")
    print(f"   Internal solver steps = {nst}")
    print(f"   Total RHS evals = {nfe}")
    print(f"   Total number of nonlinear solver iterations = {nni}")
    if args.solver == "auto":
        status, nfp, nnewt = SUNNonlinSolGetNumItersByType_Auto(NLS)
        assert status == CV_SUCCESS
        print(f"        newton={nnewt}, fixedpoint={nfp}")
    print(f"   Total number of nonlinear solver convergence failures = {ncfn}")
    if args.solver == "newton" or args.solver == "auto":
        print(f"   Total number of Jacobian evaluations = {nje}")
        print(f"   Total linear solver setups = {nsetups}")
        print(f"   Total RHS evals for setting up the linear system = {nfeLS}")

    if args.plot:
        try:
            import matplotlib.pyplot as plt
        except Exception as exc:
            print(f"\nPlot requested, but matplotlib could not be imported: {exc}")
        else:
            plt.plot(ts, y1s, label="y1")
            plt.plot(ts, y2s, label="y2")
            plt.xlabel("t")
            plt.ylabel("y")
            plt.legend()
            plt.grid(True, alpha=0.3)
            plt.title(f"Van der Pol oscillator (mu={mu:g})")
            plt.savefig(args.plot_file, dpi=150)
            print(f"\nSaved plot to: {args.plot_file}")
            if args.show:
                plt.show()
            plt.close()


def test_cvs_vdp():
    main(argv=["cvs_vdp.py"])


if __name__ == "__main__":
    main()
