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

import gc

import numpy as np
import pytest

pytest.importorskip("mpi4py")

import sundials4py.core as sun  # noqa: E402

if not all(hasattr(sun, name) for name in ("N_VMake_MPIPlusX", "N_VNew_Cuda")):
    pytest.skip("sundials4py requires MPI and CUDA support", allow_module_level=True)

cp = pytest.importorskip("cupy")
from mpi4py import MPI  # noqa: E402


@pytest.fixture
def comm():
    return MPI.COMM_WORLD


@pytest.fixture
def sunctx(comm):
    status, ctx = sun.SUNContext_Create(comm)
    assert status == sun.SUN_SUCCESS
    return ctx


@pytest.fixture(params=["device", "managed"])
def vector(request, comm, sunctx):
    n_local = 3 + comm.Get_rank()
    offset = comm.exscan(n_local) or 0
    factory = sun.N_VNew_Cuda if request.param == "device" else sun.N_VNewManaged_Cuda
    local = factory(n_local, sunctx)
    y = sun.N_VMake_MPIPlusX(comm, local, sunctx)
    sun.N_VGetCupyArray(y)[:] = offset + cp.arange(1, n_local + 1)
    cp.cuda.get_current_stream().synchronize()
    return y, local, n_local, offset


def test_global_reductions(comm, vector):
    y, local, n_local, _ = vector
    n = comm.allreduce(n_local)
    assert sun.N_VGetLength(y) == n
    assert sun.N_VGetLocalLength(y) == n_local
    assert sun.N_VGetCommunicator(y).Compare(comm) == MPI.CONGRUENT
    ones = sun.N_VClone(y)
    sun.N_VConst(1.0, ones)
    values = [
        sun.N_VDotProd(y, y),
        sun.N_VL1Norm(y),
        sun.N_VMaxNorm(y),
        sun.N_VMin(y),
        sun.N_VWrmsNorm(y, ones),
        sun.N_VMinQuotient(y, ones),
    ]
    squares = n * (n + 1) * (2 * n + 1) / 6
    np.testing.assert_allclose(
        values, [squares, n * (n + 1) / 2, n, 1, np.sqrt(squares / n), 1]
    )
    assert all(v == values for v in comm.allgather(values))


def test_local_device_views_and_clone(comm, vector):
    y, local, n_local, offset = vector
    view = sun.N_VGetCupyArray(y)
    pointer = sun.N_VGetDeviceArrayPointer(y)
    assert view.shape == (n_local,)
    assert view.data.ptr == int(pointer) == int(sun.N_VGetDeviceArrayPointer(local))
    with pytest.raises(RuntimeError, match="not supported for CUDA N_Vectors"):
        sun.N_VGetNumpyArray(y)
    z = sun.N_VClone(y)
    sun.N_VConst(2.0, z)
    sun.N_VLinearSum(1.0, y, 1.0, z, z)
    np.testing.assert_allclose(
        cp.asnumpy(sun.N_VGetCupyArray(z)), offset + np.arange(1, n_local + 1) + 2
    )
    assert sun.N_VMaxNorm(z) == pytest.approx(comm.allreduce(n_local) + 2)
    view[:] = -1.0
    cp.cuda.get_current_stream().synchronize()
    assert sun.N_VL1Norm(y) == pytest.approx(comm.allreduce(n_local))


@pytest.mark.parametrize("accessor", ["cupy", "pointer"])
def test_clone_device_view_keeps_parent_alive(vector, accessor):
    y, _, n_local, _ = vector
    clone = sun.N_VClone(y)
    sun.N_VConst(7.0, clone)
    if accessor == "cupy":
        view = sun.N_VGetCupyArray(clone)
    else:
        view = cp.from_dlpack(sun.N_VGetDeviceArrayPointer(clone))
    del clone
    gc.collect()
    np.testing.assert_allclose(cp.asnumpy(view), 7.0)
    assert view.size == n_local


def test_cuda_subcommunicator(comm):
    sub = comm.Split(comm.Get_rank() % 2)
    status, ctx = sun.SUNContext_Create(sub)
    assert status == sun.SUN_SUCCESS
    n_local = 3 + comm.Get_rank()
    local = sun.N_VNew_Cuda(n_local, ctx)
    y = sun.N_VMake_MPIPlusX(sub, local, ctx)
    sub.Free()
    # Both context and vector retain their own communicator duplicates.
    retained = sun.N_VGetCommunicator(y)
    assert sun.N_VGetLength(y) == retained.allreduce(n_local)
    sun.N_VConst(2.0, y)
    assert sun.N_VL1Norm(y) == pytest.approx(retained.allreduce(2 * n_local))


def test_mpimanyvector_mixed_serial_cuda(comm, sunctx):
    n_local = 3 + comm.Get_rank()
    host = sun.N_VNew_Serial(2, sunctx)
    device = sun.N_VNew_Cuda(n_local, sunctx)
    sun.N_VConst(1.0, host)
    sun.N_VConst(2.0, device)
    y = sun.N_VMake_MPIManyVector(comm, 2, [host, device], sunctx)
    assert sun.N_VGetLength(y) == comm.allreduce(2 + n_local)
    assert sun.N_VL1Norm(y) == pytest.approx(comm.allreduce(2 + 2 * n_local))
    assert sun.N_VDotProd(y, y) == pytest.approx(comm.allreduce(2 + 4 * n_local))
    clone = sun.N_VClone(y)
    sub = sun.N_VGetSubvector_MPIManyVector(clone, 1)
    view = sun.N_VGetCupyArray(sub)
    sun.N_VConst(3.0, clone)
    np.testing.assert_allclose(cp.asnumpy(sun.N_VGetCupyArray(sub)), 3.0)
    del clone, sub
    gc.collect()
    np.testing.assert_allclose(cp.asnumpy(view), 3.0)
