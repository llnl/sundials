#!/usr/bin/env python3
# -----------------------------------------------------------------
# Programmer(s): Cody Balos
# -----------------------------------------------------------------
# Van der Pol oscillator example using CVODE (Python)
# -----------------------------------------------------------------
import sys
import numpy as np
from sundials4py.core import *
from sundials4py.cvodes import *

class VanDerPolODE:
    def __init__(self, mu, y10, y20):
        self.mu = mu
        self.y10 = y10
        self.y20 = y20
        self.eig_log_file = "jac_eigenvalues.txt"
        # Overwrite file at start
        with open(self.eig_log_file, "w") as f:
            f.write("# t min_eig max_eig\n")

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
        ydot[1] = mu * (1.0 - y[0]**2) * y[1] - y[0]
        return 0

    def jac(self, t, yvec, fyvec, J, tmp1, tmp2, tmp3):
        y = N_VGetArrayPointer(yvec)
        mu = self.mu
        Jdata = SUNDenseMatrix_Data(J)
        Jdata[0, 0] = 0.0
        Jdata[0, 1] = 1.0
        Jdata[1, 0] = -2.0 * mu * y[0] * y[1] - 1.0
        Jdata[1, 1] = mu * (1.0 - y[0]**2)

        # Compute eigenvalues
        eigvals = np.linalg.eigvals(np.array(Jdata))
        min_eig = np.min(np.real(eigvals))
        max_eig = np.max(np.real(eigvals))
        # Log to file
        with open(self.eig_log_file, "a") as f:
            f.write(f"{t:.16e} {min_eig:.16e} {max_eig:.16e}\n")
        return 0

def main():
    mu = 10.0
    y10 = 2.0
    y20 = 0.0
    T0 = 0.0
    Tf = 300.0
    dTout = 0.01
    NEQ = 2
    Nt = int(np.ceil(Tf / dTout))
    reltol = 1e-6
    abstol = 1e-10

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

    A = SUNDenseMatrix(NEQ, NEQ, sunctx)
    LS = SUNLinSol_Dense(y, A, sunctx)

    status = CVodeSetLinearSolver(cvode.get(), LS, A)
    assert status == CV_SUCCESS

    status = CVodeSetJacFn(cvode.get(), lambda t, yvec, fyvec, J, tmp1, tmp2, tmp3, _: ode.jac(t, yvec, fyvec, J, tmp1, tmp2, tmp3))
    assert status == CV_SUCCESS

    status = CVodeSetOptions(cvode.get(), "", "", len(sys.argv), sys.argv)
    assert status == CV_SUCCESS

    yarr = N_VGetArrayPointer(y)
    print("\nVan der Pol oscillator (CVODE):")
    print(f"    initial conditions: y1 = {y10}, y2 = {y20}")
    print(f"    mu = {mu}")
    print(f"    reltol = {reltol}, abstol = {abstol}\n")
    print("        t           y1           y2")
    print("   -----------------------------------")
    print(f"  {T0:10.6f}  {yarr[0]:10.6f}  {yarr[1]:10.6f}")

    with open("cv_vdp_solution.txt", "w") as UFID:
        UFID.write("# t y1 y2\n")
        UFID.write(f" {T0:.16e} {yarr[0]:.16e} {yarr[1]:.16e}\n")
        tout = T0 + dTout
        for iout in range(Nt):
            status, tret = CVode(cvode.get(), tout, y, CV_NORMAL)
            yarr = N_VGetArrayPointer(y)
            print(f"  {tret:10.6f}  {yarr[0]:10.6f}  {yarr[1]:10.6f}")
            UFID.write(f" {tret:.16e} {yarr[0]:.16e} {yarr[1]:.16e}\n")
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
    status, nsetups = CVodeGetNumLinSolvSetups(cvode.get())
    assert status == CV_SUCCESS
    status, nni = CVodeGetNumNonlinSolvIters(cvode.get())
    assert status == CV_SUCCESS
    status, ncfn = CVodeGetNumNonlinSolvConvFails(cvode.get())
    assert status == CV_SUCCESS
    status, nje = CVodeGetNumJacEvals(cvode.get())
    assert status == CV_SUCCESS
    status, nfeLS = CVodeGetNumLinRhsEvals(cvode.get())
    assert status == CV_SUCCESS

    print("\nFinal Solver Statistics:")
    print(f"   Internal solver steps = {nst}")
    print(f"   Total RHS evals = {nfe}")
    print(f"   Total linear solver setups = {nsetups}")
    print(f"   Total number of Jacobian evaluations = {nje}")
    print(f"   Total number of Newton iterations = {nni}")
    print(f"   Total number of nonlinear solver convergence failures = {ncfn}")
    print(f"   Total RHS evals for setting up the linear system = {nfeLS}")

def test_cvs_vdp():
    main()

if __name__ == "__main__":
    main()
