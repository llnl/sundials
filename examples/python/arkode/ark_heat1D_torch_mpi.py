#!/usr/bin/env python3
# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# All rights reserved.
#
# See the top-level LICENSE and NOTICE files for details.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------
# Distributed 1D heat equation with one PyTorch CUDA N_Vector per MPI rank.
# Fixed-step explicit integration keeps halo exchanges synchronized without
# requiring an MPI N_Vector in the Python bindings.
# -----------------------------------------------------------------

import argparse
import math

import numpy as np
from mpi4py import MPI

import sundials4py.arkode as ark
import sundials4py.core as sun


def exact_semidiscrete_solution(n, k, t, dx):
    interior = np.arange(1, n - 1)
    modes = np.arange(1, n - 1)
    phi = np.sqrt(2.0 / (n - 1)) * np.sin(np.outer(interior, modes) * np.pi / (n - 1))
    lambdas = -4.0 * k / dx**2 * np.sin(0.5 * modes * np.pi / (n - 1)) ** 2
    source = np.zeros(n - 2)
    source[(n // 2) - 1] = 0.01 / dx
    solution = np.zeros(n)
    solution[1:-1] = phi @ (((np.exp(lambdas * t) - 1.0) / lambdas) * (phi.T @ source))
    return solution


def local_domain(n, rank, size):
    interior = n - 2
    count, extra = divmod(interior, size)
    local_n = count + (rank < extra)
    first = 1 + rank * count + min(rank, extra)
    return first, local_n


class DistributedHeat1D:
    def __init__(self, comm, first, local_n, n, torch, device, dx, k=0.01):
        self.comm = comm
        self.left = comm.rank - 1 if comm.rank else MPI.PROC_NULL
        self.right = comm.rank + 1 if comm.rank + 1 < comm.size else MPI.PROC_NULL
        self.source = n // 2 - first
        self.local_n = local_n
        self.torch = torch
        self.device = device
        self.coeff = k / dx**2
        self.source_strength = 0.01 / dx

    def rhs(self, t, yvec, ydotvec, user_data):
        y = sun.N_VGetTorchTensor(yvec)
        ydot = sun.N_VGetTorchTensor(ydotvec)

        # MPI transfers the two edge values through host memory. Both calls
        # use the same order on every rank, including ranks with one point.
        left = self.comm.sendrecv(
            y[-1].item(), dest=self.right, sendtag=0, source=self.left, recvtag=0
        )
        right = self.comm.sendrecv(
            y[0].item(), dest=self.left, sendtag=1, source=self.right, recvtag=1
        )

        ydot[:] = -2.0 * y
        ydot[1:] += y[:-1]
        ydot[:-1] += y[1:]
        ydot[0] += 0.0 if left is None else left
        ydot[-1] += 0.0 if right is None else right
        ydot *= self.coeff
        if 0 <= self.source < self.local_n:
            ydot[self.source] += self.source_strength
        if self.device == "cuda":
            self.torch.cuda.synchronize()
        return 0


class BatchedHeat1D:
    def __init__(self, n, fields, torch, device, dx, k=0.01):
        self.interior = n - 2
        self.fields = fields
        self.source = n // 2 - 1
        self.torch = torch
        self.device = device
        self.coeff = k / dx**2
        self.source_strength = 0.01 / dx

    def rhs(self, t, yvec, ydotvec, user_data):
        y = sun.N_VGetTorchTensor(yvec).view(self.fields, self.interior)
        ydot = sun.N_VGetTorchTensor(ydotvec).view(self.fields, self.interior)
        ydot[:] = -2.0 * y
        ydot[:, 1:] += y[:, :-1]
        ydot[:, :-1] += y[:, 1:]
        ydot *= self.coeff
        ydot[:, self.source] += self.source_strength
        if self.device == "cuda":
            self.torch.cuda.synchronize()
        return 0


def select_gpu(torch, comm):
    local_comm = comm.Split_type(MPI.COMM_TYPE_SHARED)
    local_rank = local_comm.rank
    local_comm.Free()
    visible = torch.cuda.device_count()
    if visible == 0:
        raise RuntimeError("PyTorch cannot see a CUDA device")
    # Slurm may expose one assigned GPU per task or all GPUs on the node.
    device_id = 0 if visible == 1 else local_rank
    if device_id >= visible:
        raise RuntimeError("More local MPI ranks than visible CUDA devices")
    torch.cuda.set_device(device_id)
    return device_id


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--n", type=int, default=101, help="global number of grid points"
    )
    parser.add_argument("--tf", type=float, default=1.0, help="final time")
    parser.add_argument(
        "--dx", type=float, default=None, help="fixed grid spacing for weak scaling"
    )
    parser.add_argument("--timing", action="store_true", help="report solver wall time")
    parser.add_argument(
        "--fields", type=int, default=1, help="heat problems per MPI rank"
    )
    parser.add_argument("--device", choices=("cpu", "cuda"), default="cuda")
    args = parser.parse_args()

    import torch

    comm = MPI.COMM_WORLD
    if args.fields == 1 and args.n < comm.size + 2:
        parser.error("--n must provide at least one interior point per MPI rank")
    if args.fields > 1 and args.n < 3:
        parser.error("--n must provide at least one interior point")
    if args.fields < 1:
        parser.error("--fields must be positive")
    if args.tf <= 0:
        parser.error("--tf must be positive")
    if args.dx is not None and args.dx <= 0:
        parser.error("--dx must be positive")
    if args.device == "cuda":
        if not hasattr(sun, "N_VMake_Cuda"):
            raise RuntimeError("sundials4py was not built with CUDA support")
        device_id = select_gpu(torch, comm)

    if args.fields > 1:
        local_n = args.fields * (args.n - 2)
    else:
        first, local_n = local_domain(args.n, comm.rank, comm.size)
    dx = args.dx if args.dx is not None else 1.0 / (args.n - 1)
    dtype = (
        torch.float32
        if np.dtype(sun.sunrealtype) == np.dtype(np.float32)
        else torch.float64
    )
    status, sunctx = sun.SUNContext_Create(sun.SUN_COMM_NULL)
    if status != sun.SUN_SUCCESS:
        raise RuntimeError(f"SUNContext_Create failed with status {status}")

    if args.device == "cuda":
        host_buffer = np.zeros(local_n, dtype=sun.sunrealtype)
        device_data = torch.zeros(local_n, device=f"cuda:{device_id}", dtype=dtype)
        y = sun.N_VMake_Cuda(local_n, host_buffer, device_data, sunctx)
    else:
        y = sun.N_VNew_Serial(local_n, sunctx)
    sun.N_VConst(0.0, y)
    if args.fields > 1:
        problem = BatchedHeat1D(args.n, args.fields, torch, args.device, dx)
    else:
        problem = DistributedHeat1D(
            comm, first, local_n, args.n, torch, args.device, dx
        )

    stepper = ark.ERKStepCreate(problem.rhs, 0.0, y, sunctx)
    if stepper is None:
        raise RuntimeError("ERKStepCreate failed")
    # RK4 is stable for this diffusion operator with h <= 0.4*dx^2/k.
    steps = 10 * math.ceil(args.tf / (10 * 0.4 / problem.coeff))
    step = args.tf / steps
    status = ark.ARKodeSetFixedStep(stepper.get(), step)
    if status != ark.ARK_SUCCESS:
        raise RuntimeError(f"ARKodeSetFixedStep failed with status {status}")
    status = ark.ARKodeSetMaxNumSteps(stepper.get(), steps + 10)
    if status != ark.ARK_SUCCESS:
        raise RuntimeError(f"ARKodeSetMaxNumSteps failed with status {status}")

    if comm.rank == 0:
        print(
            f"MPI + PyTorch {args.device} heat equation: "
            f"{comm.size} ranks, n={args.n} per problem, dx={dx:.3e}, "
            f"domain length={dx * (args.n - 1):.3e}, h={step:.3e}"
        )
        if args.fields > 1:
            print(f"{args.fields} fields/rank, {local_n} unknowns/GPU")
        print("        t      ||u||_rms")
    if args.timing:
        comm.Barrier()
        solve_start = MPI.Wtime()
    t = 0.0
    for output in range(1, 11):
        status, t = ark.ARKodeEvolve(
            stepper.get(), args.tf * output / 10, y, ark.ARK_NORMAL
        )
        if status != ark.ARK_SUCCESS:
            raise RuntimeError(
                f"ARKodeEvolve failed with status {status} on rank {comm.rank}"
            )
        local_y = sun.N_VGetTorchTensor(y)
        local_sum = (local_y * local_y).sum().item()
        total_n = args.n * args.fields * comm.size if args.fields > 1 else args.n
        rms = math.sqrt(comm.allreduce(local_sum, op=MPI.SUM) / total_n)
        if comm.rank == 0:
            print(f"  {t:10.6f}  {rms:10.6f}")

    if args.timing:
        solve_seconds = comm.allreduce(MPI.Wtime() - solve_start, op=MPI.MAX)
        if comm.rank == 0:
            print(f"SOLVE_SECONDS={solve_seconds:.6f}")

    local_result = sun.N_VGetTorchTensor(y).cpu().numpy().copy()
    exact = exact_semidiscrete_solution(args.n, 0.01, args.tf, dx)
    if args.fields > 1:
        local_result = local_result.reshape(args.fields, args.n - 2)
        local_error = np.max(np.abs(local_result - exact[1:-1]))
        error = comm.allreduce(local_error, op=MPI.MAX)
        if comm.rank == 0:
            print(f"Final max error vs exact semi-discrete solution = {error:.6e}")
        np.testing.assert_allclose(
            local_result,
            np.broadcast_to(exact[1:-1], local_result.shape),
            rtol=1e-4,
            atol=1e-8,
        )
        return

    parts = comm.gather(local_result, root=0)
    if comm.rank == 0:
        result = np.zeros(args.n, dtype=sun.sunrealtype)
        result[1:-1] = np.concatenate(parts)
        error = np.max(np.abs(result - exact))
        print(f"Final max error vs exact semi-discrete solution = {error:.6e}")
        np.testing.assert_allclose(result, exact, rtol=1e-4, atol=1e-8)


if __name__ == "__main__":
    main()
