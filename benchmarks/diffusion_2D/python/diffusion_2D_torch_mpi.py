#!/usr/bin/env python3
# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------
# Python/PyTorch port of the ARKODE MPI + CUDA 2D diffusion benchmark.
# -----------------------------------------------------------------

import argparse
import math

import numpy as np
from mpi4py import MPI

import sundials4py.arkode as ark
import sundials4py.core as sun


class Diffusion2D:
    def __init__(self, comm, torch, nx, rows, kx, ky, device, rhs_mode):
        self.comm = comm
        self.torch = torch
        self.nx = nx
        self.device = device
        self.rhs_mode = rhs_mode
        self.time_rhs = False
        self.rhs_seconds = 0.0
        self.rhs_calls = 0
        self.local_rows = rows + (comm.rank == 0)
        self.first_row = 0 if comm.rank == 0 else comm.rank * rows + 1
        self.south = comm.rank - 1 if comm.rank else MPI.PROC_NULL
        self.north = comm.rank + 1 if comm.rank + 1 < comm.size else MPI.PROC_NULL
        dx = 1.0 / (nx - 1)
        dy = 1.0 / rows
        self.cx = kx / (dx * dx)
        self.cy = ky / (dy * dy)
        self.cc = -2.0 * (self.cx + self.cy)
        self.inv_diag = 1.0

        dtype = (
            torch.float32 if np.dtype(sun.sunrealtype) == np.float32 else torch.float64
        )
        x = torch.arange(nx, device=device, dtype=dtype) * dx
        y = (
            torch.arange(self.local_rows, device=device, dtype=dtype) + self.first_row
        ) * dy
        sx2 = torch.sin(math.pi * x).square()[None, :]
        sy2 = torch.sin(math.pi * y).square()[:, None]
        self.shape = sx2 * sy2
        self.spatial_forcing = (
            2.0
            * math.pi**2
            * (kx * (1.0 - 2.0 * sx2) * sy2 + ky * (1.0 - 2.0 * sy2) * sx2)
        )
        self.south_send = self.south_recv = None
        self.north_send = self.north_recv = None
        if self.south != MPI.PROC_NULL:
            self.south_send = np.empty(nx, dtype=sun.sunrealtype)
            self.south_recv = np.empty(nx, dtype=sun.sunrealtype)
            self.south_send_tensor = torch.from_numpy(self.south_send)
            self.south_recv_tensor = torch.from_numpy(self.south_recv)
            self.south_recv_device = torch.empty(nx, dtype=dtype, device=device)
        if self.north != MPI.PROC_NULL:
            self.north_send = np.empty(nx, dtype=sun.sunrealtype)
            self.north_recv = np.empty(nx, dtype=sun.sunrealtype)
            self.north_send_tensor = torch.from_numpy(self.north_send)
            self.north_recv_tensor = torch.from_numpy(self.north_recv)
            self.north_recv_device = torch.empty(nx, dtype=dtype, device=device)
        if rhs_mode == "triton":
            import diffusion_2D_triton

            self.triton_rhs = diffusion_2D_triton
            self.triton_double = dtype == torch.float64
            if self.south_recv is None:
                self.south_recv_device = torch.zeros(nx, dtype=dtype, device=device)
            if self.north_recv is None:
                self.north_recv_device = torch.zeros(nx, dtype=dtype, device=device)
        self.shape[:, 0] = 0.0
        self.shape[:, -1] = 0.0
        if self.south == MPI.PROC_NULL:
            self.shape[0, :] = 0.0
        if self.north == MPI.PROC_NULL:
            self.shape[-1, :] = 0.0

    def exact(self, t):
        return self.shape * math.cos(math.pi * t) ** 2

    def tensor(self, vector):
        local = sun.N_VGetLocalVector_MPIPlusX(vector)
        return sun.N_VGetTorchTensor(local).view(self.local_rows, self.nx)

    def warmup_triton(self, y):
        ydot = self.torch.empty_like(y)
        self.south_recv_device.zero_()
        self.north_recv_device.zero_()
        words = 2 * self.triton_rhs.coefficient_words(1.0, self.triton_double)
        self.triton_rhs.launch_interior(self, y, ydot, words)
        self.triton_rhs.launch_boundary(self, y, ydot, words)
        self.torch.cuda.synchronize()

    def rhs(self, t, yvec, ydotvec, user_data):
        if self.time_rhs:
            rhs_start = MPI.Wtime()
        y = self.tensor(yvec)
        ydot = self.tensor(ydotvec)
        requests = []
        if self.south_recv is not None:
            requests.append(self.comm.Irecv(self.south_recv, source=self.south, tag=3))
        if self.north_recv is not None:
            requests.append(self.comm.Irecv(self.north_recv, source=self.north, tag=2))
        if self.south_send is not None:
            self.south_send_tensor.copy_(y[0, :])
            requests.append(self.comm.Isend(self.south_send, dest=self.south, tag=2))
        if self.north_send is not None:
            self.north_send_tensor.copy_(y[-1, :])
            requests.append(self.comm.Isend(self.north_send, dest=self.north, tag=3))

        angle = math.pi * t
        if self.rhs_mode == "inplace":
            self.torch.mul(y, self.cc, out=ydot)
            ydot[:, 1:].add_(y[:, :-1], alpha=self.cx)
            ydot[:, :-1].add_(y[:, 1:], alpha=self.cx)
            ydot[1:, :].add_(y[:-1, :], alpha=self.cy)
            ydot[:-1, :].add_(y[1:, :], alpha=self.cy)
            ydot.add_(
                self.shape,
                alpha=-2.0 * math.pi * math.sin(angle) * math.cos(angle),
            )
            ydot.add_(self.spatial_forcing, alpha=-math.cos(angle) ** 2)
        elif self.rhs_mode == "triton":
            words = (
                *self.triton_rhs.coefficient_words(
                    -2.0 * math.pi * math.sin(angle) * math.cos(angle),
                    self.triton_double,
                ),
                *self.triton_rhs.coefficient_words(
                    -math.cos(angle) ** 2, self.triton_double
                ),
            )
            self.triton_rhs.launch_interior(self, y, ydot, words)
        else:
            ydot[:] = 0.0
            ydot[:, 1:-1] = self.cx * (y[:, :-2] - 2.0 * y[:, 1:-1] + y[:, 2:])
            ydot[1:-1, :] += self.cy * (y[:-2, :] - 2.0 * y[1:-1, :] + y[2:, :])
            ydot += (
                -2.0 * math.pi * math.sin(angle) * math.cos(angle) * self.shape
                - math.cos(angle) ** 2 * self.spatial_forcing
            )
        if requests:
            MPI.Request.Waitall(requests)
        if self.south_recv is not None:
            self.south_recv_device.copy_(self.south_recv_tensor)
            if self.rhs_mode == "inplace":
                ydot[0, :].add_(self.south_recv_device, alpha=self.cy)
            elif self.rhs_mode == "tensor":
                ydot[0, :] += self.cy * (
                    self.south_recv_device - 2.0 * y[0, :] + y[1, :]
                )
        if self.north_recv is not None:
            self.north_recv_device.copy_(self.north_recv_tensor)
            if self.rhs_mode == "inplace":
                ydot[-1, :].add_(self.north_recv_device, alpha=self.cy)
            elif self.rhs_mode == "tensor":
                ydot[-1, :] += self.cy * (
                    y[-2, :] - 2.0 * y[-1, :] + self.north_recv_device
                )
        if self.rhs_mode == "triton":
            self.triton_rhs.launch_boundary(self, y, ydot, words)
        else:
            ydot[:, 0] = 0.0
            ydot[:, -1] = 0.0
            if self.south == MPI.PROC_NULL:
                ydot[0, :] = 0.0
            if self.north == MPI.PROC_NULL:
                ydot[-1, :] = 0.0
        # The Triton kernels, halo copies, and default CUDA N_Vector policy all
        # use the default stream, so SUNDIALS sees the completed RHS in order.
        if self.device.startswith("cuda") and (
            self.rhs_mode != "triton" or self.time_rhs
        ):
            self.torch.cuda.synchronize()
        if self.time_rhs:
            self.rhs_seconds += MPI.Wtime() - rhs_start
            self.rhs_calls += 1
        return 0

    def psetup(self, t, y, fy, jok, gamma, user_data):
        self.inv_diag = 1.0 / (1.0 + 2.0 * gamma * (self.cx + self.cy))
        return 0, True

    def psolve(self, t, y, fy, r, z, gamma, delta, lr, user_data):
        sun.N_VScale(self.inv_diag, r, z)
        return 0


def select_gpu(torch, comm):
    local = comm.Split_type(MPI.COMM_TYPE_SHARED)
    local_rank = local.rank
    local.Free()
    visible = torch.cuda.device_count()
    if visible == 0:
        raise RuntimeError("PyTorch cannot see a CUDA GPU")
    device_id = 0 if visible == 1 else local_rank
    if device_id >= visible:
        raise RuntimeError("More local ranks than visible CUDA GPUs")
    torch.cuda.set_device(device_id)
    return f"cuda:{device_id}"


def check(status, name):
    if status != ark.ARK_SUCCESS:
        raise RuntimeError(f"{name} failed: {status}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nx", type=int, default=32)
    parser.add_argument("--rows-per-rank", type=int, default=31)
    parser.add_argument("--kx", type=float, default=1.0)
    parser.add_argument("--ky", type=float, default=1.0)
    parser.add_argument("--tf", type=float, default=1.0)
    parser.add_argument("--nout", type=int, default=20)
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cuda")
    parser.add_argument(
        "--rhs-mode", choices=("tensor", "inplace", "triton"), default="tensor"
    )
    parser.add_argument(
        "--time-rhs",
        action="store_true",
        help="time RHS callbacks with CUDA synchronization",
    )
    args = parser.parse_args()
    if args.nx < 3 or args.rows_per_rank < 2:
        parser.error("--nx must be at least 3 and --rows-per-rank at least 2")
    if args.kx <= 0 or args.ky <= 0 or args.tf <= 0 or args.nout < 1:
        parser.error("--kx, --ky, --tf and --nout must be positive")
    if args.rhs_mode == "triton" and args.device != "cuda":
        parser.error("--rhs-mode triton requires --device cuda")

    import torch

    comm = MPI.COMM_WORLD
    if not hasattr(sun, "N_VMake_MPIPlusX"):
        raise RuntimeError("sundials4py was not built with MPIPlusX support")
    if args.device == "cuda" and not hasattr(sun, "N_VMake_Cuda"):
        raise RuntimeError("sundials4py was not built with CUDA support")
    device = select_gpu(torch, comm) if args.device == "cuda" else "cpu"
    problem = Diffusion2D(
        comm,
        torch,
        args.nx,
        args.rows_per_rank,
        args.kx,
        args.ky,
        device,
        args.rhs_mode,
    )
    local_n = args.nx * problem.local_rows
    dtype = torch.float32 if np.dtype(sun.sunrealtype) == np.float32 else torch.float64
    status, ctx = sun.SUNContext_Create(sun.SUN_COMM_NULL)
    if status != sun.SUN_SUCCESS:
        raise RuntimeError(f"SUNContext_Create failed: {status}")
    if args.device == "cuda":
        host = np.zeros(local_n, dtype=sun.sunrealtype)
        data = torch.zeros(local_n, dtype=dtype, device=device)
        local_vector = sun.N_VMake_Cuda(local_n, host, data, ctx)
    else:
        local_vector = sun.N_VNew_Serial(local_n, ctx)
    y = sun.N_VMake_MPIPlusX(comm.py2f(), local_vector, ctx)
    if y is None:
        raise RuntimeError("N_VMake_MPIPlusX failed")
    problem.tensor(y)[:] = problem.exact(0.0)
    if args.device == "cuda":
        torch.cuda.synchronize()
    if args.rhs_mode == "triton":
        problem.warmup_triton(problem.tensor(y))

    stepper = ark.ARKStepCreate(None, problem.rhs, 0.0, y, ctx)
    if stepper is None:
        raise RuntimeError("ARKStepCreate failed")
    mem = stepper.get()
    check(ark.ARKodeSStolerances(mem, 1e-5, 1e-10), "ARKodeSStolerances")
    linear_solver = sun.SUNLinSol_PCG(y, sun.SUN_PREC_RIGHT, 20, ctx)
    if linear_solver is None:
        raise RuntimeError("SUNLinSol_PCG failed")
    check(ark.ARKodeSetLinearSolver(mem, linear_solver, None), "ARKodeSetLinearSolver")
    check(
        ark.ARKodeSetPreconditioner(mem, problem.psetup, problem.psolve),
        "ARKodeSetPreconditioner",
    )
    check(ark.ARKodeSetLSetupFrequency(mem, 0), "ARKodeSetLSetupFrequency")
    check(ark.ARKodeSetEpsLin(mem, 0.0), "ARKodeSetEpsLin")
    check(ark.ARKodeSetOrder(mem, 3), "ARKodeSetOrder")
    check(ark.ARKodeSetLinear(mem, 0), "ARKodeSetLinear")
    check(
        ark.ARKodeSetAdaptControllerByName(mem, "I"), "ARKodeSetAdaptControllerByName"
    )
    check(ark.ARKodeSetMaxNumSteps(mem, 0), "ARKodeSetMaxNumSteps")
    check(ark.ARKodeSetStopTime(mem, args.tf), "ARKodeSetStopTime")

    if comm.rank == 0:
        print(
            f"2D diffusion ARKStep: {comm.size} ranks, nx={args.nx}, "
            f"ny={args.rows_per_rank * comm.size + 1}, "
            f"kx={args.kx}, ky={args.ky}, tf={args.tf}, nout={args.nout}, "
            f"rhs_mode={args.rhs_mode}"
        )
    comm.Barrier()
    problem.time_rhs = args.time_rhs
    solve_seconds = 0.0
    max_error = 0.0
    final_error = 0.0
    for output in range(1, args.nout + 1):
        tout = args.tf * output / args.nout
        start = MPI.Wtime()
        status, t = ark.ARKodeEvolve(mem, tout, y, ark.ARK_NORMAL)
        if args.device == "cuda":
            torch.cuda.synchronize()
        solve_seconds += MPI.Wtime() - start
        if status == ark.ARK_TSTOP_RETURN:
            if output != args.nout or not math.isclose(t, args.tf, abs_tol=1e-12):
                raise RuntimeError(f"ARKodeEvolve reached stop time early: t={t}")
        else:
            check(status, "ARKodeEvolve")
        error = (problem.tensor(y) - problem.exact(t)).abs().max().item()
        final_error = comm.allreduce(error, op=MPI.MAX)
        max_error = max(max_error, final_error)
    if args.time_rhs:
        rank_times = comm.allgather(
            (solve_seconds, problem.rhs_seconds, problem.rhs_calls)
        )
        timing_rank = max(range(comm.size), key=lambda rank: rank_times[rank][0])
        solve_seconds, rhs_seconds, rhs_calls = rank_times[timing_rank]
    else:
        solve_seconds = comm.allreduce(solve_seconds, op=MPI.MAX)
    stats = {
        "STEPS": ark.ARKodeGetNumSteps,
        "STEP_ATTEMPTS": ark.ARKodeGetNumStepAttempts,
        "ERROR_TEST_FAILS": ark.ARKodeGetNumErrTestFails,
        "NLS_ITERS": ark.ARKodeGetNumNonlinSolvIters,
        "NLS_FAILS": ark.ARKodeGetNumNonlinSolvConvFails,
        "LS_ITERS": ark.ARKodeGetNumLinIters,
        "LS_FAILS": ark.ARKodeGetNumLinConvFails,
        "PREC_SETUPS": ark.ARKodeGetNumPrecEvals,
        "PREC_SOLVES": ark.ARKodeGetNumPrecSolves,
    }
    counters = {}
    for name, getter in stats.items():
        status, value = getter(mem)
        check(status, name)
        counters[name] = comm.allreduce(value, op=MPI.MAX)
    if comm.rank == 0:
        print(f"SOLVE_SECONDS={solve_seconds:.6f}")
        if args.time_rhs:
            print(f"TIMING_RANK={timing_rank}")
            print(f"RHS_CALLS={rhs_calls}")
            print(f"RHS_SECONDS={rhs_seconds:.6f}")
            print(f"NON_RHS_SECONDS={solve_seconds - rhs_seconds:.6f}")
            print(f"RHS_PERCENT={100.0 * rhs_seconds / solve_seconds:.2f}")
        print(f"FINAL_ERROR={final_error:.6e}")
        print(f"MAX_ERROR={max_error:.6e}")
        for name, value in counters.items():
            print(f"{name}={value}")
    if not math.isfinite(max_error) or max_error > 0.01:
        raise RuntimeError(f"Solution error exceeds 0.01: {max_error}")


if __name__ == "__main__":
    main()
