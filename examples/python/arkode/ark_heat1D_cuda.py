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
# This example is a copy of examples/arkode/C_serial/ark_heat1D.c,
# but ported to Python to use sundials4py as well as CUDA through
# either CuPy or PyTorch.
#
# The following test simulates a simple 1D heat equation,
#    u_t = k*u_xx + f
# for t in [0, 10], x in [0, 1], with initial conditions
#    u(0,x) =  0
# Dirichlet boundary conditions, i.e.
#    u_t(t,0) = u_t(t,1) = 0,
# and a point-source heating term,
#    f = 0.01 for x=0.5.
#
# The spatial derivatives are computed using second-order
# centered differences, with the data distributed over N points
# on a uniform spatial grid.
#
# The final solution is checked against the exact solution of this
# semi-discrete ODE system.  With zero Dirichlet boundaries, the
# interior finite-difference Laplacian is diagonalized by the discrete
# sine basis.  For zero initial data and a constant point source b, the
# interior solution is
#
#    u(t) = Phi * diag((exp(lambda_m*t) - 1) / lambda_m) * Phi^T * b,
#
# where Phi contains the orthonormal sine modes and lambda_m are the
# corresponding finite-difference Laplacian eigenvalues.  This validates
# the time integration error without introducing error from a continuous
# PDE approximation.
#
# This program solves the problem with either an ERK or DIRK
# method.  For the DIRK method, we use a Newton iteration with
# the SUNLinSol_PCG linear solver, and a user-supplied Jacobian-vector
# product routine.
#
# 100 outputs are printed at equal intervals, and run statistics
# are printed at the end.
# -----------------------------------------------------------------

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


class Heat1DCudaProblem:
    def __init__(self, backend, n=101, k=0.01):
        self.backend = backend
        self.n = n
        self.k = k
        self.dx = 1.0 / (n - 1)
        self.isource = n // 2

    def device_array(self, nvec):
        return self.backend.from_dlpack(sun.N_VGetDeviceArray_Cuda(nvec))

    def set_init_cond(self, yvec):
        y = self.device_array(yvec)
        y[:] = 0.0

    def f(self, t, yvec, ydotvec, user_data):
        y = self.device_array(yvec)
        ydot = self.device_array(ydotvec)

        ydot[:] = 0.0
        c1 = self.k / self.dx / self.dx
        c2 = -2.0 * self.k / self.dx / self.dx
        ydot[1:-1] = c1 * y[:-2] + c2 * y[1:-1] + c1 * y[2:]
        ydot[0] = 0.0
        ydot[-1] = 0.0
        ydot[self.isource] += 0.01 / self.dx
        self.backend.synchronize()
        return 0


class CupyBackend:
    name = "cupy"

    def __init__(self):
        import cupy

        self.xp = cupy

    def zeros(self, n):
        return self.xp.zeros(n, dtype=sun.sunrealtype)

    def from_dlpack(self, obj):
        return self.xp.from_dlpack(obj)

    def to_numpy(self, array):
        return self.xp.asnumpy(array)

    def synchronize(self):
        self.xp.cuda.Device().synchronize()


class TorchBackend:
    name = "torch"

    def __init__(self):
        import torch

        if not torch.cuda.is_available():
            raise RuntimeError("PyTorch CUDA is not available")
        self.torch = torch
        self.dtype = torch.float32 if sun.sunrealtype == np.float32 else torch.float64

    def zeros(self, n):
        return self.torch.zeros(n, device="cuda", dtype=self.dtype)

    def from_dlpack(self, obj):
        return self.torch.utils.dlpack.from_dlpack(obj)

    def to_numpy(self, array):
        return array.cpu().numpy()

    def synchronize(self):
        self.torch.cuda.synchronize()


def solve_heat1d(backend):
    n = 101
    k = 0.01
    tf = 1.0
    nt = 10
    reltol = 1e-6
    abstol = 1e-10

    status, sunctx = sun.SUNContext_Create(sun.SUN_COMM_NULL)
    assert status == sun.SUN_SUCCESS

    host_data = np.zeros(n, dtype=sun.sunrealtype)
    device_data = backend.zeros(n)
    y = sun.N_VMake_Cuda(n, host_data, device_data, sunctx)

    problem = Heat1DCudaProblem(backend, n=n, k=k)
    problem.set_init_cond(y)

    stepper = ark.ARKStepCreate(problem.f, None, 0.0, y, sunctx)
    assert stepper is not None

    status = ark.ARKodeSStolerances(stepper.get(), reltol, abstol)
    assert status == ark.ARK_SUCCESS

    status = ark.ARKodeSetMaxNumSteps(stepper.get(), 100000)
    assert status == ark.ARK_SUCCESS

    t = 0.0
    tout = tf / nt

    print(f"\n{backend.name} backend")
    print("        t      ||u||_rms")
    print("   -------------------------")
    print(f"  {t:10.6f}  {0.0:10.6f}")

    for _ in range(nt):
        status, t = ark.ARKodeEvolve(stepper.get(), tout, y, ark.ARK_NORMAL)
        if status != ark.ARK_SUCCESS:
            raise RuntimeError(f"ARKodeEvolve failed with status {status}")

        sun.N_VCopyFromDevice_Cuda(y)
        rms = np.sqrt(np.dot(host_data, host_data) / n)
        print(f"  {t:10.6f}  {rms:10.6f}")
        tout = min(tout + tf / nt, tf)

    uexact = exact_semidiscrete_solution(n, k, tf)
    max_error = np.max(np.abs(host_data - uexact))
    print(f"\nFinal max error vs exact semi-discrete solution = {max_error:.6e}")
    np.testing.assert_allclose(host_data, uexact, rtol=1e-4, atol=1e-8)

    return host_data.copy(), backend.to_numpy(device_data)


def main():
    if not hasattr(sun, "N_VMake_Cuda"):
        raise RuntimeError("sundials4py was not built with CUDA support")

    results = {}
    for backend_type in (CupyBackend, TorchBackend):
        try:
            backend = backend_type()
        except ImportError:
            continue
        except RuntimeError as err:
            print(f"Skipping {backend_type.name}: {err}")
            continue

        host_data, device_data = solve_heat1d(backend)
        np.testing.assert_allclose(host_data, device_data)
        results[backend.name] = host_data

    if not results:
        raise RuntimeError("Install CuPy or PyTorch with CUDA support to run this example")

    if {"cupy", "torch"} <= results.keys():
        np.testing.assert_allclose(results["cupy"], results["torch"], rtol=1e-6, atol=1e-10)


if __name__ == "__main__":
    main()
