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
# either CuPy, PyTorch, or JAX.
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


class Heat1DCudaProblem:
    def __init__(self, backend, n=101, k=0.01):
        self.backend = backend
        self.n = n
        self.k = k
        self.dx = 1.0 / (n - 1)
        self.isource = n // 2

    def f(self, t, yvec, ydotvec, user_data):
        y = self.backend.get_array(yvec)

        self.backend.heat_rhs(y, ydotvec, self.k, self.dx, self.isource)
        return 0

    def jtv(self, vvec, Jvvec, t, yvec, fyvec, user_data, tmpvec):
        V = self.backend.get_array(vvec)

        self.backend.heat_jtv(V, Jvvec, self.k, self.dx)
        return 0


# CuPy and PyTorch expose mutable CUDA arrays, so the SUNDIALS output vectors
# can be filled in place by the callback functions below.
class MutableCudaArrayBackend:
    def set_zero(self, nvec):
        array = self.get_array(nvec)
        array[:] = 0.0
        self.synchronize()

    def heat_rhs(self, y, ydotvec, k, dx, isource):
        ydot = self.get_array(ydotvec)
        ydot[:] = 0.0
        c1 = k / dx / dx
        c2 = -2.0 * k / dx / dx
        ydot[1:-1] = c1 * y[:-2] + c2 * y[1:-1] + c1 * y[2:]
        ydot[0] = 0.0
        ydot[-1] = 0.0
        ydot[isource] += 0.01 / dx
        self.synchronize()

    def heat_jtv(self, v, Jvvec, k, dx):
        Jv = self.get_array(Jvvec)
        Jv[:] = 0.0
        c1 = k / dx / dx
        c2 = -2.0 * k / dx / dx
        Jv[1:-1] = c1 * v[:-2] + c2 * v[1:-1] + c1 * v[2:]
        Jv[0] = 0.0
        Jv[-1] = 0.0
        self.synchronize()


class CupyBackend(MutableCudaArrayBackend):
    name = "cupy"
    get_array = staticmethod(sun.N_VGetCupyArray)

    def __init__(self):
        import cupy

        self.xp = cupy
        self.synchronize = cupy.cuda.Device().synchronize


class TorchBackend(MutableCudaArrayBackend):
    name = "torch"
    get_array = staticmethod(sun.N_VGetTorchTensor)

    def __init__(self):
        import torch

        if not torch.cuda.is_available():
            raise RuntimeError("PyTorch CUDA is not available")
        self.torch = torch
        self.dtype = torch.float32 if sun.sunrealtype == np.float32 else torch.float64
        self.synchronize = torch.cuda.synchronize


# JAX arrays are immutable. The JAX callbacks compute a new device array and
# then replace the N_Vector device pointer with N_VSetDeviceArrayPointer_Cuda.
class JaxBackend:
    name = "jax"
    get_array = staticmethod(sun.N_VGetJaxArray)

    def __init__(self):
        import jax

        if np.dtype(sun.sunrealtype) == np.dtype(np.float64):
            jax.config.update("jax_enable_x64", True)

        import jax.numpy as jnp

        devices = [
            device for device in jax.devices() if device.platform in ("cuda", "gpu")
        ]
        if not devices:
            raise RuntimeError("JAX CUDA/GPU backend is not available")

        self.jax = jax
        self.xp = jnp
        self.device = devices[0]
        self.dtype = jnp.float32 if sun.sunrealtype == np.float32 else jnp.float64

    def set_zero(self, nvec):
        array = self.xp.zeros(sun.N_VGetLength(nvec), dtype=self.dtype)
        array = self.jax.device_put(array, self.device).block_until_ready()
        sun.N_VSetDeviceArrayPointer_Cuda(array, nvec)

    def heat_rhs(self, y, ydotvec, k, dx, isource):
        c1 = k / dx / dx
        c2 = -2.0 * k / dx / dx
        result = self.xp.zeros_like(y)
        result = result.at[1:-1].set(c1 * y[:-2] + c2 * y[1:-1] + c1 * y[2:])
        result = result.at[isource].add(0.01 / dx)
        sun.N_VSetDeviceArrayPointer_Cuda(result.block_until_ready(), ydotvec)

    def heat_jtv(self, v, Jvvec, k, dx):
        c1 = k / dx / dx
        c2 = -2.0 * k / dx / dx
        result = self.xp.zeros_like(v)
        result = result.at[1:-1].set(c1 * v[:-2] + c2 * v[1:-1] + c1 * v[2:])
        sun.N_VSetDeviceArrayPointer_Cuda(result.block_until_ready(), Jvvec)


def solve_heat1d(backend, device_data, n=101):
    k = 0.01
    tf = 1.0
    nt = 10
    reltol = 1e-6
    abstol = 1e-10

    status, sunctx = sun.SUNContext_Create(sun.SUN_COMM_NULL)
    assert status == sun.SUN_SUCCESS

    host_data = np.zeros(n, dtype=sun.sunrealtype)
    y = sun.N_VMake_Cuda(n, host_data, device_data, sunctx)

    problem = Heat1DCudaProblem(backend, n=n, k=k)
    backend.set_zero(y)

    stepper = ark.ARKStepCreate(None, problem.f, 0.0, y, sunctx)
    assert stepper is not None

    status = ark.ARKodeSStolerances(stepper.get(), reltol, abstol)
    assert status == ark.ARK_SUCCESS

    status = ark.ARKodeSetMaxNumSteps(stepper.get(), 100000)
    assert status == ark.ARK_SUCCESS

    # PCG linear solver with no preconditioning, with up to n iterations
    LS = sun.SUNLinSol_PCG(y, sun.SUN_PREC_NONE, n, sunctx)

    status = ark.ARKodeSetLinearSolver(stepper.get(), LS, None)
    assert status == ark.ARK_SUCCESS

    status = ark.ARKodeSetJacTimes(stepper.get(), None, problem.jtv)
    assert status == ark.ARK_SUCCESS

    status = ark.ARKodeSetLinear(stepper.get(), 0)
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

    # Print statistics
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--n", type=int, default=101, help="number of spatial grid points"
    )
    args = parser.parse_args()
    if args.n < 3:
        parser.error("--n must be at least 3")

    if not hasattr(sun, "N_VMake_Cuda"):
        raise RuntimeError("sundials4py was not built with CUDA support")

    results = {}
    for backend_type in (CupyBackend, TorchBackend, JaxBackend):
        try:
            backend = backend_type()
        except ImportError:
            continue
        except RuntimeError as err:
            print(f"Skipping {backend_type.name}: {err}")
            continue

        n = args.n
        if isinstance(backend, CupyBackend):
            device_data = backend.xp.zeros(n, dtype=sun.sunrealtype)
            host_data, y = solve_heat1d(backend, device_data, n)
            device_result = backend.xp.asnumpy(backend.get_array(y))
        elif isinstance(backend, TorchBackend):
            device_data = backend.torch.zeros(n, device="cuda", dtype=backend.dtype)
            host_data, y = solve_heat1d(backend, device_data, n)
            device_result = backend.get_array(y).cpu().numpy()
        else:
            device_data = backend.xp.zeros(n, dtype=backend.dtype)
            device_data = backend.jax.device_put(device_data, backend.device)
            host_data, y = solve_heat1d(backend, device_data, n)
            device_result = np.asarray(backend.get_array(y))

        np.testing.assert_allclose(host_data, device_result)
        results[backend.name] = host_data

    if not results:
        raise RuntimeError(
            "Install CuPy, PyTorch, or JAX with CUDA support to run this example"
        )

    reference_name = next(iter(results))
    for name, result in results.items():
        np.testing.assert_allclose(
            results[reference_name], result, rtol=1e-6, atol=1e-10
        )


if __name__ == "__main__":
    main()
