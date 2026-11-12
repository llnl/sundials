#!/usr/bin/env python3
# -----------------------------------------------------------------
# Programmer(s): Cody J. Balos @ LLNL
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
# This example solves the 1D heat equation
#   u_t = k u_xx,  u(0,t) = u(1,t) = 0,
# distributed over MPI ranks with an MPIPlusX N_Vector over a serial local
# vector, mpi4py halo exchange on NumPy views, and ARKStep with SPGMR
# (classical Gram-Schmidt, which needs fewer global reductions).
#
# Requires a sundials4py build with SUNDIALS_ENABLE_MPI=ON and mpi4py.
# Run with:
#   mpiexec -n 4 python -m mpi4py ark_heat1D_mpi.py
# --fields F solves F independent heat equations on the same spatial grid
# (default F=1). Each grid point stores F values, so the global vector has nx*F
# entries and each rank stores (local grid rows, F). Halo exchanges send F values
# per boundary row. Increasing F increases vector work and halo message size
# without refining the grid or changing the diffusion operator's stiffness.
# By default, every field starts with u(x,0) = sin(pi*x).
# With --multimode, field j starts with
#   u_j(x,0) = a_j*sin(pi*x) + b_j*sin(3*pi*x),
#   a_j = 1 + 0.25*cos(2*pi*j/fields),
#   b_j = 0.25 + 0.1*sin(2*pi*j/fields),  j = 0,...,fields-1.
# The sin(3*pi*x) component decays faster under diffusion. --multimode
# changes only the initial condition; fields remain independent and the PDE,
# grid, and MPI decomposition are unchanged.
# Use --solver sts for stabilized explicit diffusion integration without Krylov
# solves. The JAX and Torch versions share this numerical and timing driver.
# -----------------------------------------------------------------

import argparse
import json
import time

import numpy as np

try:
    from mpi4py import MPI
except ImportError:
    MPI = None

import sundials4py.arkode as ark
import sundials4py.core as sun

K = 0.5
TF = 0.1
RTOL = 1.0e-6
ATOL = 1.0e-10
ERROR_LIMIT = 2.0e-5


def mpi_available():
    return MPI is not None and hasattr(sun, "N_VMake_MPIPlusX")


def check_status(status):
    if status != 0:
        raise RuntimeError(f"SUNDIALS returned status {status}")


class Heat1D:
    """Independent heat fields, stored as (local grid rows, fields)."""

    backend = "numpy"
    kernel = "numpy"

    def __init__(self, comm, nx, fields, eager=False):
        self.comm = comm
        self.nx = nx
        self.fields = fields
        self.multimode = False
        rank, size = comm.Get_rank(), comm.Get_size()
        base, rem = divmod(nx, size)
        self.n_local = base + (rank < rem)
        self.offset = rank * base + min(rank, rem)
        self.shape = (self.n_local, fields)
        self.dx = 1.0 / (nx + 1)
        self.coefficient = K / self.dx**2
        self.left = rank - 1 if rank > 0 else MPI.PROC_NULL
        self.right = rank + 1 if rank < size - 1 else MPI.PROC_NULL
        self.halo = np.zeros((2, fields), dtype=sun.sunrealtype)
        self.ext = np.empty((self.n_local + 2, fields), dtype=sun.sunrealtype)

    def exact(self, t):
        # Exact solution of the discretized Laplacian, not the continuum PDE.
        # Each initial sine mode decays as exp(eigenvalue*t); t=0 initializes y.
        x = (self.offset + np.arange(1, self.n_local + 1)) * self.dx
        q = np.arange(self.fields) / self.fields
        a = (
            1.0 + 0.25 * np.cos(2 * np.pi * q)
            if self.multimode
            else np.ones(self.fields)
        )
        b = (
            0.25 + 0.1 * np.sin(2 * np.pi * q)
            if self.multimode
            else np.zeros(self.fields)
        )
        solution = np.zeros(self.shape, dtype=sun.sunrealtype)
        for mode, amplitude in ((1, a), (3, b)):
            eigenvalue = -4 * self.coefficient * np.sin(mode * np.pi * self.dx / 2) ** 2
            solution += (
                np.exp(eigenvalue * t) * np.sin(mode * np.pi * x[:, None]) * amplitude
            )
        return solution

    def exchange(self, boundary):
        """Exchange completed host rows; PROC_NULL leaves zero boundary values."""
        self.halo.fill(0)
        self.comm.Sendrecv(
            np.ascontiguousarray(boundary[1]),
            dest=self.right,
            recvbuf=self.halo[0],
            source=self.left,
            sendtag=0,
            recvtag=0,
        )
        self.comm.Sendrecv(
            np.ascontiguousarray(boundary[0]),
            dest=self.left,
            recvbuf=self.halo[1],
            source=self.right,
            sendtag=1,
            recvtag=1,
        )
        return self.halo

    def new_local(self, ctx):
        return sun.N_VNew_Serial(self.n_local * self.fields, ctx)

    def initialize(self, vector, values):
        sun.N_VGetNumpyArray(vector)[:] = values.ravel()

    def host(self, vector):
        return np.array(sun.N_VGetNumpyArray(vector)).reshape(self.shape)

    def synchronize(self, vector=None):
        pass

    def apply(self, vector, output):
        u = sun.N_VGetNumpyArray(vector).reshape(self.shape)
        halo = self.exchange(u[[0, -1]])
        self.ext[1:-1] = u
        self.ext[0] = halo[0]
        self.ext[-1] = halo[1]
        f = sun.N_VGetNumpyArray(output).reshape(self.shape)
        f[:] = self.coefficient * (self.ext[:-2] - 2 * u + self.ext[2:])

    def f(self, t, yvec, fvec, user_data):
        self.apply(yvec, fvec)
        return 0

    def jtv(self, vvec, jvvec, t, yvec, fyvec, user_data, tmpvec):
        self.apply(vvec, jvvec)
        return 0

    def metadata(self):
        return {
            "backend": self.backend,
            "kernel": self.kernel,
            "numpy_version": np.__version__,
        }

    def dom_eig(self, t, yvec, fvec, user_data, tmp1, tmp2, tmp3):
        return 0, -4 * self.coefficient, 0.0


def solve_once(problem, ctx, solver):
    setup_start = time.perf_counter()
    local = problem.new_local(ctx)
    y = sun.N_VMake_MPIPlusX(problem.comm, local, ctx)
    if sun.N_VGetLength(y) != problem.nx * problem.fields:
        raise RuntimeError("Distributed vector has the wrong global length")
    problem.initialize(y, problem.exact(0))
    if solver == "arkstep":
        integrator = ark.ARKStepCreate(None, problem.f, 0.0, y, ctx)
        ls = sun.SUNLinSol_SPGMR(y, 0, 20, ctx)
        check_status(sun.SUNLinSol_SPGMRSetGSType(ls, sun.SUN_CLASSICAL_GS))
        check_status(ark.ARKodeSetLinearSolver(integrator.get(), ls, None))
        check_status(ark.ARKodeSetJacTimes(integrator.get(), None, problem.jtv))
        check_status(ark.ARKodeSetLinear(integrator.get(), 0))
    else:
        integrator = ark.LSRKStepCreateSTS(problem.f, 0.0, y, ctx)
        check_status(ark.LSRKStepSetSTSMethod(integrator.get(), ark.ARKODE_LSRK_RKC_2))
        check_status(ark.LSRKStepSetDomEigFn(integrator.get(), problem.dom_eig))
        check_status(ark.LSRKStepSetDomEigFrequency(integrator.get(), 25))
        check_status(ark.LSRKStepSetDomEigSafetyFactor(integrator.get(), 1.01))
        check_status(ark.LSRKStepSetMaxNumStages(integrator.get(), 1000))
        check_status(ark.ARKodeSetAdaptControllerByName(integrator.get(), "I"))
    check_status(ark.ARKodeSStolerances(integrator.get(), RTOL, ATOL))
    check_status(ark.ARKodeSetMaxNumSteps(integrator.get(), 10000))
    problem.synchronize(y)
    setup_seconds = problem.comm.allreduce(
        time.perf_counter() - setup_start, op=MPI.MAX
    )

    problem.comm.Barrier()
    start = MPI.Wtime()
    status, t = ark.ARKodeEvolve(integrator.get(), TF, y, ark.ARK_NORMAL)
    problem.synchronize(y)
    elapsed = MPI.Wtime() - start
    solve_seconds = problem.comm.allreduce(elapsed, op=MPI.MAX)
    check_status(status)
    if not np.isclose(t, TF):
        raise RuntimeError("Solver did not reach the requested final time")
    error = problem.host(y) - problem.exact(TF)
    max_error = problem.comm.allreduce(float(np.max(np.abs(error))), op=MPI.MAX)
    if not np.isfinite(max_error) or max_error >= ERROR_LIMIT:
        raise RuntimeError(f"Analytic solution error {max_error} exceeds {ERROR_LIMIT}")
    counts = {}
    getters = [
        ("steps", ark.ARKodeGetNumSteps),
        ("step_attempts", ark.ARKodeGetNumStepAttempts),
        ("error_test_failures", ark.ARKodeGetNumErrTestFails),
    ]
    implicit_getters = [
        ("linear_iters", ark.ARKodeGetNumLinIters),
        ("jv_evals", ark.ARKodeGetNumJtimesEvals),
        ("linear_failures", ark.ARKodeGetNumLinConvFails),
        ("nonlinear_failures", ark.ARKodeGetNumNonlinSolvConvFails),
    ]
    if solver == "arkstep":
        getters.extend(implicit_getters)
    else:
        counts.update({name: 0 for name, _ in implicit_getters})
    for name, getter in getters:
        flag, count = getter(integrator.get())
        check_status(flag)
        counts[name] = count
    flag, counts["rhs_evals"] = ark.ARKodeGetNumRhsEvals(
        integrator.get(), 1 if solver == "arkstep" else 0
    )
    check_status(flag)
    if counts["steps"] <= 0 or counts["rhs_evals"] <= 0:
        raise RuntimeError("No integration work completed")
    if not all(c == counts for c in problem.comm.allgather(counts)):
        raise RuntimeError("Solver work differs between ranks")
    return dict(
        status=status,
        final_time=t,
        solve_seconds=solve_seconds,
        setup_seconds=setup_seconds,
        setup_and_solve_seconds=setup_seconds + solve_seconds,
        max_error=max_error,
        error_limit=ERROR_LIMIT,
        **counts,
    )


def main(argv=None, problem_type=Heat1D):
    parser = argparse.ArgumentParser(
        description="Distributed independent 1D heat equations"
    )
    parser.add_argument("--nx", type=int, default=201, help="global interior grid rows")
    parser.add_argument(
        "--fields",
        type=int,
        default=1,
        help="number of independent heat equations on the same grid (nx*fields unknowns)",
    )
    parser.add_argument(
        "--json", action="store_true", help="print one machine-readable result"
    )
    parser.add_argument(
        "--eager", action="store_true", help="disable Torch stencil compilation"
    )
    parser.add_argument(
        "--solver",
        choices=("arkstep", "sts"),
        default="arkstep",
        help="implicit ARKStep or explicit stabilized LSRKStep",
    )
    parser.add_argument(
        "--multimode",
        action="store_true",
        help="initialize each field with sin(pi*x) and sin(3*pi*x), with field-dependent amplitudes",
    )
    args = parser.parse_args(argv)
    if not mpi_available():
        parser.error("This example needs mpi4py and an MPI-enabled sundials4py build")
    comm = MPI.COMM_WORLD
    if args.nx < 3 or args.fields < 1 or comm.Get_size() > args.nx:
        parser.error("Require --nx >= 3, --fields >= 1, and at most --nx ranks")
    if args.eager and problem_type.backend != "torch":
        parser.error("--eager applies only to the Torch example")
    start = time.perf_counter()
    problem = problem_type(comm, args.nx, args.fields, args.eager)
    problem.multimode = args.multimode
    status, ctx = sun.SUNContext_Create(comm)
    check_status(status)
    preparation_seconds = comm.allreduce(time.perf_counter() - start, op=MPI.MAX)
    start = time.perf_counter()
    solve_once(
        problem, ctx, args.solver
    )  # Warm the actual callback, compiler and solver paths.
    warmup_seconds = comm.allreduce(time.perf_counter() - start, op=MPI.MAX)
    result = solve_once(problem, ctx, args.solver)
    result.update(
        nx=args.nx,
        fields=args.fields,
        ranks=comm.Get_size(),
        k=K,
        tf=TF,
        rtol=RTOL,
        atol=ATOL,
        dtype=np.dtype(sun.sunrealtype).name,
        solver=(
            "ARKStep implicit, linear, SPGMR(20), classical GS, analytic Jv"
            if args.solver == "arkstep"
            else "LSRKStep RKC2, analytic spectral bound, I controller"
        ),
        multimode=args.multimode,
        preparation_seconds=preparation_seconds,
        warmup_seconds=warmup_seconds,
        mpi_version=MPI.Get_library_version().strip("\x00\n"),
        **problem.metadata(),
    )
    if comm.Get_rank() == 0:
        if args.json:
            print(json.dumps(result), flush=True)
        else:
            print(
                f"1D heat equation: {args.nx} rows, {args.fields} fields, {comm.Get_size()} ranks"
            )
            print(f"   backend = {result['backend']}, kernel = {result['kernel']}")
            print(
                f"   steps = {result['steps']}, linear iterations = {result['linear_iters']}"
            )
            print(f"   max error = {result['max_error']:.3e}")
            print(f"   warmed solve = {result['solve_seconds']:.6f} s")
            print(
                f"   setup = {result['setup_seconds']:.6f} s, warmup = {warmup_seconds:.6f} s"
            )
    return result


def test_ark_heat1d_mpi():
    if mpi_available():
        main([])


def test_ark_heat1d_mpi_sts():
    if mpi_available():
        main(["--nx", "7", "--fields", "3", "--solver", "sts", "--multimode"])


if __name__ == "__main__":
    main()
