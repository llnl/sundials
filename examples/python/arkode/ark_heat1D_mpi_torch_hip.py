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
# MPIPlusX(HIP) version of ark_heat1D_mpi.py using mutable Torch tensors.
# By default each independent field starts with sin(pi*x). --multimode adds
# sin(3*pi*x) and varies both amplitudes by field, as defined in ark_heat1D_mpi.py.
# The sin(3*pi*x) component decays faster. Only the initial condition changes.
# --fields F batches F independent heat equations on the same grid, giving nx*F
# unknowns. Each rank stores F values per local grid row and exchanges F values
# per boundary row. Larger F increases vector work without refining the grid.
# Requires MPI and HIP bindings, mpi4py, and ROCm PyTorch.
# Halo rows pass through host buffers; HIP-aware MPI is not required.
# Mask ROCR_VISIBLE_DEVICES per rank to choose GPUs on multi-GPU systems.
# The stencil is compiled by default. Pass --eager to use eager Torch instead.
#   mpiexec -n 2 python -m mpi4py ark_heat1D_mpi_torch_hip.py --fields 256
# -----------------------------------------------------------------

import numpy as np

import sundials4py.core as sun

if __package__:
    from .ark_heat1D_mpi import Heat1D, main
else:
    from ark_heat1D_mpi import Heat1D, main


class TorchHeat1D(Heat1D):
    backend = "torch"

    def __init__(self, comm, nx, fields, eager=False):
        super().__init__(comm, nx, fields, eager)
        if not hasattr(sun, "N_VNew_Hip"):
            raise RuntimeError("sundials4py was built without HIP support")
        import torch

        if not torch.cuda.is_available() or torch.version.hip is None:
            raise RuntimeError("PyTorch ROCm is not available")
        torch.set_num_threads(1)
        torch.cuda.set_device(0)
        self.torch = torch
        self.dtype = (
            torch.float64
            if np.dtype(sun.sunrealtype) == np.dtype(np.float64)
            else torch.float32
        )

        def stencil(u, halo, output):
            ext = torch.cat((halo[:1], u, halo[1:]), dim=0)
            output.copy_(self.coefficient * (ext[:-2] - 2 * u + ext[2:]))

        self.stencil = stencil if eager else torch.compile(stencil, fullgraph=True)
        self.kernel = "torch eager" if eager else "torch.compile"

    def new_local(self, ctx):
        return sun.N_VNew_Hip(self.n_local * self.fields, ctx)

    def initialize(self, vector, values):
        sun.N_VGetTorchTensor_Hip(vector).copy_(
            self.torch.as_tensor(values.ravel(), device="cuda", dtype=self.dtype)
        )
        self.synchronize()

    def host(self, vector):
        return (
            sun.N_VGetTorchTensor_Hip(vector).cpu().numpy().copy().reshape(self.shape)
        )

    def synchronize(self, vector=None):
        self.torch.cuda.synchronize()

    def apply(self, vector, output):
        with self.torch.no_grad():
            u = sun.N_VGetTorchTensor_Hip(vector).reshape(self.shape)
            boundary = self.torch.stack((u[0], u[-1])).cpu().numpy()
            halo = self.torch.as_tensor(
                self.exchange(boundary), device="cuda", dtype=self.dtype
            )
            f = sun.N_VGetTorchTensor_Hip(output).reshape(self.shape)
            self.stencil(u, halo, f)
            # SUNDIALS does not run on Torch's stream. Finish writes before return.
            self.synchronize()

    def metadata(self):
        return dict(
            **super().metadata(),
            framework_version=self.torch.__version__,
            hip_version=self.torch.version.hip,
            device="cuda:0",
            device_kind=self.torch.cuda.get_device_name(0),
        )


def test_ark_heat1d_mpi_torch():
    import pytest

    if not hasattr(sun, "N_VMake_MPIPlusX") or not hasattr(sun, "N_VNew_Hip"):
        pytest.skip("MPI and HIP bindings are required")
    pytest.importorskip("mpi4py")
    torch = pytest.importorskip("torch")
    if not torch.cuda.is_available() or torch.version.hip is None:
        pytest.skip("Torch ROCm device is not available")
    main(["--nx", "7", "--fields", "3", "--solver", "sts", "--multimode"], TorchHeat1D)


if __name__ == "__main__":
    main(problem_type=TorchHeat1D)
