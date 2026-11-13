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
# MPIPlusX(HIP) version of ark_heat1D_mpi.py using immutable JAX arrays.
# By default each independent field starts with sin(pi*x). --multimode adds
# sin(3*pi*x) and varies both amplitudes by field, as defined in ark_heat1D_mpi.py.
# The sin(3*pi*x) component decays faster. Only the initial condition changes.
# --fields F batches F independent heat equations on the same grid, giving nx*F
# unknowns. Each rank stores F values per local grid row and exchanges F values
# per boundary row. Larger F increases vector work without refining the grid.
# Requires MPI and HIP bindings, mpi4py, and JAX with a ROCm backend.
# Halo rows pass through host buffers; HIP-aware MPI is not required.
# For shared GPUs, set XLA_PYTHON_CLIENT_PREALLOCATE=false before launching.
# Mask ROCR_VISIBLE_DEVICES per rank to choose GPUs on multi-GPU systems.
#   mpiexec -n 2 python -m mpi4py ark_heat1D_mpi_jax_hip.py --fields 256
# -----------------------------------------------------------------

import numpy as np

import sundials4py.core as sun

if __package__:
    from .ark_heat1D_mpi import Heat1D, main
else:
    from ark_heat1D_mpi import Heat1D, main


class JaxHeat1D(Heat1D):
    backend = "jax"
    kernel = "jax.jit"

    def __init__(self, comm, nx, fields, eager=False):
        super().__init__(comm, nx, fields, eager)
        if not hasattr(sun, "N_VNew_Hip"):
            raise RuntimeError("sundials4py was built without HIP support")
        import jax

        if np.dtype(sun.sunrealtype) == np.dtype(np.float64):
            jax.config.update("jax_enable_x64", True)
        import jax.numpy as jnp

        self.jax = jax
        self.jnp = jnp
        devices = [
            device
            for device in jax.devices()
            if device.platform == "rocm"
            or (
                device.platform == "gpu"
                and "rocm" in device.client.platform_version.lower()
            )
        ]
        if not devices:
            raise RuntimeError("JAX ROCm backend is not available")
        self.device = devices[0]
        # Select the same visible GPU in the HIP runtime before vector allocation.
        jax.device_put(
            jnp.zeros(1, dtype=sun.sunrealtype), self.device
        ).block_until_ready()

        def stencil(u, halo):
            ext = jnp.concatenate((halo[:1], u, halo[1:]), axis=0)
            return (self.coefficient * (ext[:-2] - 2 * u + ext[2:])).reshape(-1)

        self.stencil = jax.jit(stencil)
        self.boundary = jax.jit(lambda u: jnp.stack((u[0], u[-1])))

    def new_local(self, ctx):
        return sun.N_VNew_Hip(self.n_local * self.fields, ctx)

    def initialize(self, vector, values):
        sun.N_VSetJaxArray_Hip(self.jax.device_put(values.ravel(), self.device), vector)

    def host(self, vector):
        return np.array(self.jax.device_get(sun.N_VGetJaxArray_Hip(vector))).reshape(
            self.shape
        )

    def synchronize(self, vector=None):
        self.jax.effects_barrier()
        if vector is not None:
            # This accessor waits on SUNDIALS' HIP stream, including interpolation.
            sun.N_VGetJaxArray_Hip(vector).block_until_ready()

    def apply(self, vector, output):
        u = sun.N_VGetJaxArray_Hip(vector).reshape(self.shape)
        boundary = np.asarray(self.jax.device_get(self.boundary(u)))
        halo = self.jax.device_put(self.exchange(boundary), self.device)
        result = self.stencil(u, halo)
        # JAX arrays are immutable. The setter waits and copies into solver storage.
        sun.N_VSetJaxArray_Hip(result, output)

    def metadata(self):
        return dict(
            **super().metadata(),
            framework_version=self.jax.__version__,
            device=str(self.device),
            device_kind=self.device.device_kind,
        )


def test_ark_heat1d_mpi_jax():
    import pytest

    from sundials4py import core

    if not hasattr(core, "N_VMake_MPIPlusX") or not hasattr(core, "N_VNew_Hip"):
        pytest.skip("MPI and HIP bindings are required")
    pytest.importorskip("mpi4py")
    jax = pytest.importorskip("jax")
    if not any(
        device.platform == "rocm"
        or (
            device.platform == "gpu"
            and "rocm" in device.client.platform_version.lower()
        )
        for device in jax.devices()
    ):
        pytest.skip("JAX ROCm device is not available")
    main(["--nx", "7", "--fields", "3", "--solver", "sts", "--multimode"], JaxHeat1D)


if __name__ == "__main__":
    main(problem_type=JaxHeat1D)
