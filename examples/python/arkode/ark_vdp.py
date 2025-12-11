#!/usr/bin/env python3
# -----------------------------------------------------------------
# Programmer(s): Cody Balos
# -----------------------------------------------------------------
# Van der Pol oscillator example using ARKODE (Python)
# -----------------------------------------------------------------
import numpy as np
from sundials4py.core import *
from sundials4py.arkode import *


class VanDerPolODE:
    def __init__(self, mu, y10, y20, sunctx):
        self.sunctx = sunctx
        self.ark = None
        self.mu = mu
        self.y10 = y10
        self.y20 = y20
        self.eig_log_file = "jac_eigenvalues.txt"
        # Overwrite file at start
        with open(self.eig_log_file, "w") as f:
            f.write("# t stiff_test\n")

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

        J = SUNDenseMatrix(2, 2, self.sunctx)
        self.jac(t, yvec, ydotvec, J, None, None, None)

        return 0

    def jac(self, t, yvec, fyvec, J, tmp1, tmp2, tmp3):
        y = N_VGetArrayPointer(yvec)
        mu = self.mu
        Jdata = SUNDenseMatrix_Data(J)
        Jdata[0, 0] = 0.0
        Jdata[0, 1] = 1.0
        Jdata[1, 0] = -2.0 * mu * y[0] * y[1] - 1.0
        Jdata[1, 1] = mu * (1.0 - y[0] ** 2)

        # Compute eigenvalues
        stab = 3.5068  # Tsit5
        eigvals = np.linalg.eigvals(np.array(Jdata))
        max_eig = np.max(np.real(eigvals))
        status, h = ARKodeGetCurrentStep(self.ark.get())
        # Log to file
        with open(self.eig_log_file, "a") as f:
            f.write(f"{t:.16e} {np.abs(h*max_eig)/stab:.16e}\n")
        return 0


def main():
    mu = 10.0
    y10 = 2.0
    y20 = 0.0
    T0 = 0.0
    Tf = 150.0
    dTout = 0.001
    NEQ = 2
    Nt = int(np.ceil(Tf / dTout))
    reltol = 1e-4
    abstol = 1e-12

    status, sunctx = SUNContext_Create(SUN_COMM_NULL)
    y = N_VNew_Serial(NEQ, sunctx)

    ode = VanDerPolODE(mu, y10, y20, sunctx)
    ode.set_init_cond(y)

    implicit = False

    if implicit:
        ark = ARKStepCreate(
            None, lambda t, yvec, ydotvec, _: ode.f(t, yvec, ydotvec), T0, y, sunctx
        )
    else:
        ark = ARKStepCreate(
            lambda t, yvec, ydotvec, _: ode.f(t, yvec, ydotvec), None, T0, y, sunctx
        )

        ARKStepSetTableName(ark.get(), "", "ARKODE_TSITOURAS_7_4_5")

    ode.ark = ark

    status = ARKodeSStolerances(ark.get(), reltol, abstol)
    assert status == ARK_SUCCESS

    status = ARKodeSetInterpolantType(ark.get(), ARK_INTERP_LAGRANGE)
    assert status == ARK_SUCCESS

    if implicit:
        status = ARKodeSetDeduceImplicitRhs(ark.get(), 1)
        assert status == ARK_SUCCESS

        A = SUNDenseMatrix(NEQ, NEQ, sunctx)
        LS = SUNLinSol_Dense(y, A, sunctx)

        status = ARKodeSetLinearSolver(ark.get(), LS, A)
        assert status == ARK_SUCCESS

        status = ARKodeSetJacFn(
            ark.get(),
            lambda t, yvec, fyvec, J, tmp1, tmp2, tmp3, _: ode.jac(
                t, yvec, fyvec, J, tmp1, tmp2, tmp3
            ),
        )
        assert status == ARK_SUCCESS

        status = ARKodeSetAutonomous(ark.get(), 1)
        assert status == ARK_SUCCESS

    yarr = N_VGetArrayPointer(y)
    print("\nVan der Pol oscillator (ARKODE):")
    print(f"    initial conditions: y1 = {y10}, y2 = {y20}")
    print(f"    mu = {mu}")
    print(f"    reltol = {reltol}, abstol = {abstol}\n")
    print("        t           y1           y2")
    print("   -----------------------------------")
    print(f"  {T0:10.6f}  {yarr[0]:10.6f}  {yarr[1]:10.6f}")

    with open("ark_vdp_solution.txt", "w") as UFID:
        UFID.write("# t y1 y2\n")
        UFID.write(f" {T0:.16e} {yarr[0]:.16e} {yarr[1]:.16e}\n")
        tout = T0 + dTout
        for iout in range(Nt):
            status, tret = ARKodeEvolve(ark.get(), tout, y, ARK_NORMAL)
            yarr = N_VGetArrayPointer(y)
            print(f"  {tret:10.6f}  {yarr[0]:10.6f}  {yarr[1]:10.6f}")
            UFID.write(f" {tret:.16e} {yarr[0]:.16e} {yarr[1]:.16e}\n")
            if status == ARK_SUCCESS:
                tout += dTout
                tout = min(tout, Tf)
            else:
                print("Solver failure, stopping integration")
                break
        print("   -----------------------------------")

    status, nst = ARKodeGetNumSteps(ark.get())
    assert status == ARK_SUCCESS
    status, nst_a = ARKodeGetNumStepAttempts(ark.get())
    assert status == ARK_SUCCESS
    status, nfe = ARKodeGetNumRhsEvals(ark.get(), 0)
    assert status == ARK_SUCCESS
    status, nfi = ARKodeGetNumRhsEvals(ark.get(), 1)
    assert status == ARK_SUCCESS
    status, netf = ARKodeGetNumErrTestFails(ark.get())
    assert status == ARK_SUCCESS
    status, ncfn = ARKodeGetNumStepSolveFails(ark.get())
    assert status == ARK_SUCCESS
    if implicit:
        status, nsetups = ARKodeGetNumLinSolvSetups(ark.get())
        assert status == ARK_SUCCESS
        status, nni = ARKodeGetNumNonlinSolvIters(ark.get())
        assert status == ARK_SUCCESS
        status, nnf = ARKodeGetNumNonlinSolvConvFails(ark.get())
        assert status == ARK_SUCCESS
        status, nje = ARKodeGetNumJacEvals(ark.get())
        assert status == ARK_SUCCESS
        status, nfeLS = ARKodeGetNumLinRhsEvals(ark.get())
        assert status == ARK_SUCCESS

    print("\nFinal Solver Statistics:")
    print(f"   Internal solver steps = {nst} (attempted = {nst_a})")
    print(f"   Total RHS evals:  Fe = {nfe},  Fi = {nfi}")
    print(f"   Total number of error test failures = {netf}")
    if implicit:
        print(f"   Total linear solver setups = {nsetups}")
        print(f"   Total RHS evals for setting up the linear system = {nfeLS}")
        print(f"   Total number of Jacobian evaluations = {nje}")
        print(f"   Total number of Newton iterations = {nni}")
        print(f"   Total number of nonlinear solver convergence failures = {nnf}")
        print(f"   Total number of failed steps from solver failure = {ncfn}")


def test_ark_vdp():
    main()


if __name__ == "__main__":
    main()
