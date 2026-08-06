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
# 1D heat equation example using sundials4py, CUDA N_Vectors, and CuPy.
# -----------------------------------------------------------------

import argparse

import numpy as np
import sundials4py.arkode as ark
import sundials4py.core as sun


def exact_semidiscrete_solution(n, k, t):
    dx = 1.0 / (n - 1)
    i = np.arange(1, n - 1)
    m = np.arange(1, n - 1)

    phi = np.sqrt(2.0 / (n - 1)) * np.sin(np.outer(i, m) * np.pi / (n - 1))
    lambdas = -4.0 * k / dx**2 * np.sin(0.5 * m * np.pi / (n - 1)) ** 2

    source = np.zeros(n - 2)
    source[(n // 2) - 1] = 0.01 / dx

    u = np.zeros(n)
    u[1:-1] = phi @ (((np.exp(lambdas * t) - 1.0) / lambdas) * (phi.T @ source))
    return u


def solve_heat1d(name, y, host_data, problem, sunctx, n=101, k=0.01):
    tf = 1.0
    nt = 10
    reltol = 1e-6
    abstol = 1e-10

    problem.set_init_cond(y)

    stepper = ark.ARKStepCreate(None, problem.f, 0.0, y, sunctx)
    assert stepper is not None

    status = ark.ARKodeSStolerances(stepper.get(), reltol, abstol)
    assert status == ark.ARK_SUCCESS

    status = ark.ARKodeSetMaxNumSteps(stepper.get(), 100000)
    assert status == ark.ARK_SUCCESS

    # PCG linear solver with no preconditioning, with up to n iterations.
    LS = sun.SUNLinSol_PCG(y, sun.SUN_PREC_NONE, n, sunctx)

    status = ark.ARKodeSetLinearSolver(stepper.get(), LS, None)
    assert status == ark.ARK_SUCCESS

    status = ark.ARKodeSetJacTimes(stepper.get(), None, problem.jtv)
    assert status == ark.ARK_SUCCESS

    status = ark.ARKodeSetLinear(stepper.get(), 0)
    assert status == ark.ARK_SUCCESS

    t = 0.0
    tout = tf / nt

    print(f"\n{name}")
    print("        t      ||u||_rms")
    print("   -------------------------")
    print(f"  {t:10.6f}  {0.0:10.6f}")

    for _ in range(nt):
        status, t = ark.ARKodeEvolve(stepper.get(), tout, y, ark.ARK_NORMAL)
        if status != ark.ARK_SUCCESS:
            raise RuntimeError(f"ARKodeEvolve failed with status {status}")

        host_data = sun.N_VGetNumpyArray(y, device="cpu", copy_from="device")
        rms = np.sqrt(np.dot(host_data, host_data) / n)
        print(f"  {t:10.6f}  {rms:10.6f}")
        tout = min(tout + tf / nt, tf)

    uexact = exact_semidiscrete_solution(n, k, tf)
    max_error = np.max(np.abs(host_data - uexact))
    print(f"\nFinal max error vs exact semi-discrete solution = {max_error:.6e}")
    np.testing.assert_allclose(host_data, uexact, rtol=1e-4, atol=1e-8)

    status, nst = ark.ARKodeGetNumSteps(stepper.get())
    assert status == ark.ARK_SUCCESS
    status, nst_a = ark.ARKodeGetNumStepAttempts(stepper.get())
    assert status == ark.ARK_SUCCESS
    status, nfe = ark.ARKodeGetNumRhsEvals(stepper.get(), 0)
    assert status == ark.ARK_SUCCESS
    status, nfi = ark.ARKodeGetNumRhsEvals(stepper.get(), 1)
    assert status == ark.ARK_SUCCESS
    status, nsetups = ark.ARKodeGetNumLinSolvSetups(stepper.get())
    assert status == ark.ARK_SUCCESS
    status, nli = ark.ARKodeGetNumLinIters(stepper.get())
    assert status == ark.ARK_SUCCESS
    status, nJv = ark.ARKodeGetNumJtimesEvals(stepper.get())
    assert status == ark.ARK_SUCCESS
    status, nlcf = ark.ARKodeGetNumLinConvFails(stepper.get())
    assert status == ark.ARK_SUCCESS
    status, nni = ark.ARKodeGetNumNonlinSolvIters(stepper.get())
    assert status == ark.ARK_SUCCESS
    status, ncfn = ark.ARKodeGetNumNonlinSolvConvFails(stepper.get())
    assert status == ark.ARK_SUCCESS
    status, netf = ark.ARKodeGetNumErrTestFails(stepper.get())
    assert status == ark.ARK_SUCCESS

    print("\nFinal Solver Statistics:")
    print(f"   Internal solver steps = {nst} (attempted = {nst_a})")
    print(f"   Total RHS evals:  Fe = {nfe},  Fi = {nfi}")
    print(f"   Total linear solver setups = {nsetups}")
    print(f"   Total linear iterations = {nli}")
    print(f"   Total number of Jacobian-vector products = {nJv}")
    print(f"   Total number of linear solver convergence failures = {nlcf}")
    print(f"   Total number of Newton iterations = {nni}")
    print(f"   Total number of nonlinear solver convergence failures = {ncfn}")
    print(f"   Total number of error test failures = {netf}")

    return host_data.copy(), y


# CuPy exposes mutable CUDA arrays, so the SUNDIALS output vectors can be
# filled in place by the callback functions below.
class CupyHeat1DProblem:
    def __init__(self, cupy, n=101, k=0.01):
        self.cupy = cupy
        self.n = n
        self.k = k
        self.dx = 1.0 / (n - 1)
        self.isource = n // 2

    def set_init_cond(self, yvec):
        y = sun.N_VGetCupyArray(yvec)
        y[:] = 0.0
        self.cupy.cuda.Device().synchronize()

    def f(self, t, yvec, ydotvec, user_data):
        y = sun.N_VGetCupyArray(yvec)
        ydot = sun.N_VGetCupyArray(ydotvec)

        ydot[:] = 0.0
        c1 = self.k / self.dx / self.dx
        c2 = -2.0 * self.k / self.dx / self.dx
        ydot[1:-1] = c1 * y[:-2] + c2 * y[1:-1] + c1 * y[2:]
        ydot[0] = 0.0
        ydot[-1] = 0.0
        ydot[self.isource] += 0.01 / self.dx
        self.cupy.cuda.Device().synchronize()
        return 0

    def jtv(self, vvec, Jvvec, t, yvec, fyvec, user_data, tmpvec):
        v = sun.N_VGetCupyArray(vvec)
        Jv = sun.N_VGetCupyArray(Jvvec)

        Jv[:] = 0.0
        c1 = self.k / self.dx / self.dx
        c2 = -2.0 * self.k / self.dx / self.dx
        Jv[1:-1] = c1 * v[:-2] + c2 * v[1:-1] + c1 * v[2:]
        Jv[0] = 0.0
        Jv[-1] = 0.0
        self.cupy.cuda.Device().synchronize()
        return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=101, help="number of spatial grid points")
    args = parser.parse_args()
    if args.n < 3:
        parser.error("--n must be at least 3")

    if not hasattr(sun, "N_VMake_Cuda"):
        raise RuntimeError("sundials4py was not built with CUDA support")

    import cupy

    status, sunctx = sun.SUNContext_Create(sun.SUN_COMM_NULL)
    assert status == sun.SUN_SUCCESS

    host_data = np.zeros(args.n, dtype=sun.sunrealtype)
    device_data = cupy.zeros(args.n, dtype=sun.sunrealtype)
    y = sun.N_VMake_Cuda(args.n, host_data, device_data, sunctx)

    problem = CupyHeat1DProblem(cupy, n=args.n)
    host_result, y = solve_heat1d("cupy cuda backend", y, host_data, problem, sunctx, n=args.n)
    device_result = cupy.asnumpy(sun.N_VGetCupyArray(y))
    np.testing.assert_allclose(host_result, device_result)


if __name__ == "__main__":
    main()
