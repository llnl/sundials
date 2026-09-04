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
import weakref

import numpy as np
import pytest
from numpy.testing import assert_allclose
from fixtures import *
from sundials4py.core import *


def test_nvector_array_helpers_are_always_available():
    assert "N_VGetNumpyArray" in globals()
    assert "N_VGetJaxArray" in globals()
    assert "N_VGetCupyArray" in globals()
    assert "N_VGetTorchTensor" in globals()
    assert "N_VSetJaxArray" in globals()
    assert "N_VSetHostArrayPointer_Cuda" not in globals()
    assert "N_VSetDeviceArrayPointer_Cuda" not in globals()
    assert "N_VGetHostArrayPointer_Cuda" not in globals()
    assert "N_VGetDeviceArrayPointer_Cuda" not in globals()


def test_create_manyvector(sunctx):
    y = N_VNew_Serial(5, sunctx)
    z = N_VNew_Serial(5, sunctx)

    N_VConst(1.0, y)
    N_VConst(2.0, z)

    yz = N_VNew_ManyVector(2, [y, z], sunctx)

    yarr = N_VGetArrayPointer(N_VGetSubvector_ManyVector(yz, 0))
    assert_allclose(N_VGetArrayPointer(y), 1.0)

    zarr = N_VGetArrayPointer(N_VGetSubvector_ManyVector(yz, 1))
    assert_allclose(N_VGetArrayPointer(z), 2.0)

    N_VConst(3.0, yz)
    assert_allclose(3.0, yarr)
    assert_allclose(3.0, zarr)


@pytest.mark.parametrize("vector_type", ["serial"])
def test_create_nvector(vector_type, sunctx):
    if vector_type == "serial":
        nvec = N_VNew_Serial(5, sunctx)
    else:
        raise ValueError("Unknown vector type")
    assert nvec is not None

    arr = N_VGetArrayPointer(nvec)
    assert arr.shape[0] == 5

    arr[:] = np.array([5.0, 4.0, 3.0, 2.0, 1.0], dtype=sunrealtype)
    assert_allclose(N_VGetArrayPointer(nvec), [5.0, 4.0, 3.0, 2.0, 1.0])

    N_VConst(2.0, nvec)
    assert_allclose(arr, 2.0)


@pytest.mark.parametrize("vector_type", ["serial"])
def test_make_nvector(vector_type, sunctx):
    arr = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=sunrealtype)
    if vector_type == "serial":
        nvec = N_VMake_Serial(5, arr, sunctx)
    else:
        raise ValueError("Unknown vector type")
    assert nvec is not None

    assert_allclose(N_VGetArrayPointer(nvec), arr)

    N_VConst(2.0, nvec)
    assert_allclose(arr, 2.0)

    arr[:] = np.array([5.0, 4.0, 3.0, 2.0, 1.0], dtype=sunrealtype)
    assert_allclose(N_VGetArrayPointer(nvec), [5.0, 4.0, 3.0, 2.0, 1.0])


def test_make_nvector_serial_retains_array(sunctx):
    array = np.arange(5, dtype=sunrealtype)
    array_ref = weakref.ref(array)
    nvec = N_VMake_Serial(5, array, sunctx)

    del array
    gc.collect()
    assert array_ref() is not None

    N_VConst(3.0, nvec)
    assert_allclose(array_ref(), 3.0)

    del nvec
    gc.collect()
    assert array_ref() is None


def test_nvector_array_helpers_serial(sunctx):
    nvec = N_VNew_Serial(5, sunctx)
    N_VConst(1.0, nvec)

    arr = N_VGetNumpyArray(nvec)
    assert_allclose(arr, 1.0)

    arr[:] = np.arange(5, dtype=sunrealtype)
    assert_allclose(N_VGetArrayPointer(nvec), np.arange(5, dtype=sunrealtype))
    assert_allclose(N_VGetNumpyArray(nvec), arr)

    with pytest.raises(RuntimeError):
        N_VGetCupyArray(nvec)
    with pytest.raises(TypeError):
        N_VGetNumpyArray(nvec, device="cuda")
    with pytest.raises(TypeError):
        N_VGetTorchTensor(nvec, device="cuda")


def test_nvector_array_helpers_serial_jax(sunctx):
    jax = pytest.importorskip("jax")

    nvec = N_VNew_Serial(5, sunctx)
    N_VConst(3.0, nvec)

    array = N_VGetJaxArray(nvec)
    assert array.device.platform == "cpu"
    assert_allclose(np.asarray(array), 3.0)

    N_VConst(4.0, nvec)
    array = N_VGetJaxArray(nvec)
    assert_allclose(np.asarray(array), 4.0)

    with pytest.raises(TypeError):
        N_VGetJaxArray(nvec, device="cpu")


def test_set_nvector_serial_jax_array(sunctx):
    jax = pytest.importorskip("jax")
    if np.dtype(sunrealtype) == np.dtype(np.float64):
        jax.config.update("jax_enable_x64", True)

    import jax.numpy as jnp

    nvec = N_VNew_Serial(5, sunctx)
    original_pointer = N_VGetArrayPointer(nvec).ctypes.data
    array = jax.device_put(
        jnp.arange(5, dtype=_jax_dtype(jnp)), jax.devices("cpu")[0]
    ).block_until_ready()

    N_VSetJaxArray(array, nvec)
    assert N_VGetArrayPointer(nvec).ctypes.data == original_pointer
    assert_allclose(N_VGetArrayPointer(nvec), np.arange(5, dtype=sunrealtype))

    N_VConst(6.0, nvec)
    assert_allclose(np.asarray(array), np.arange(5, dtype=sunrealtype))

    with pytest.raises(RuntimeError, match="copy=False requires a jax.ref.Ref"):
        N_VSetJaxArray(array, nvec, copy=False)
    with pytest.raises(RuntimeError, match="one-dimensional"):
        N_VSetJaxArray(jnp.zeros((1, 5), dtype=_jax_dtype(jnp)), nvec)
    with pytest.raises(RuntimeError, match="vector length"):
        N_VSetJaxArray(jnp.zeros(4, dtype=_jax_dtype(jnp)), nvec)
    with pytest.raises(RuntimeError, match="floating-point dtype"):
        N_VSetJaxArray(jnp.arange(5, dtype=jnp.int32), nvec)
    with pytest.raises(RuntimeError, match="SUNDIALS precision"):
        N_VSetJaxArray(jnp.zeros(5, dtype=jnp.float16), nvec)


def test_set_nvector_serial_jax_ref(sunctx):
    jax = pytest.importorskip("jax")
    if np.dtype(sunrealtype) == np.dtype(np.float64):
        jax.config.update("jax_enable_x64", True)

    import jax.numpy as jnp

    array = jax.device_put(
        jnp.arange(5, dtype=_jax_dtype(jnp)), jax.devices("cpu")[0]
    ).block_until_ready()
    ref = jax.ref.new_ref(array)
    ref[...].block_until_ready()
    nvec = N_VNew_Serial(5, sunctx)

    ref_pointer = ref.unsafe_buffer_pointer()
    ref_weak = weakref.ref(ref)
    N_VSetJaxArray(ref, nvec)
    assert N_VGetArrayPointer(nvec).ctypes.data == ref_pointer

    with pytest.raises(RuntimeError, match="copy=True requires a jax.Array"):
        N_VSetJaxArray(ref, nvec, copy=True)

    del array, ref
    gc.collect()
    assert ref_weak() is not None

    N_VConst(5.0, nvec)
    assert_allclose(np.asarray(ref_weak()[...]), 5.0)

    del nvec
    gc.collect()
    assert ref_weak() is None


def test_nvector_array_helpers_serial_torch(sunctx):
    torch = pytest.importorskip("torch")

    nvec = N_VNew_Serial(5, sunctx)
    N_VConst(2.0, nvec)

    tensor = N_VGetTorchTensor(nvec)
    assert_allclose(tensor.numpy(), 2.0)

    tensor[:] = torch.arange(5, dtype=tensor.dtype)
    assert_allclose(N_VGetArrayPointer(nvec), np.arange(5, dtype=sunrealtype))


@pytest.mark.skipif("N_VNew_Cuda" not in globals(), reason="CUDA bindings are not enabled")
def test_create_nvector_cuda(sunctx):
    nvec = _cuda_nvector_or_fail(lambda: N_VNew_Cuda(5, sunctx))

    assert N_VGetLength(nvec) == 5
    assert int(N_VGetDeviceArrayPointer(nvec)) != 0

    N_VConst(2.0, nvec)


def _cuda_nvector_or_fail(factory):
    try:
        nvec = factory()
    except RuntimeError as err:
        pytest.fail(f"CUDA vector allocation failed: {err}")

    if nvec is None:
        pytest.fail("CUDA vector allocation failed")
    return nvec


def _check_cuda_nvector_const(nvec, value):
    assert N_VGetLength(nvec) == 5
    assert int(N_VGetDeviceArrayPointer(nvec)) != 0

    N_VConst(value, nvec)


@pytest.mark.skipif("N_VNewManaged_Cuda" not in globals(), reason="CUDA bindings are not enabled")
def test_create_nvector_cuda_managed(sunctx):
    nvec = _cuda_nvector_or_fail(lambda: N_VNewManaged_Cuda(5, sunctx))
    _check_cuda_nvector_const(nvec, 3.0)


@pytest.mark.skipif(
    "N_VNewWithMemHelp_Cuda" not in globals(), reason="CUDA bindings are not enabled"
)
@pytest.mark.parametrize("use_managed_mem", [False, True])
def test_create_nvector_cuda_with_memhelp(use_managed_mem, sunctx):
    mem_helper = SUNMemoryHelper_Cuda(sunctx)
    if mem_helper is None:
        pytest.fail("CUDA memory helper creation failed")

    nvec = _cuda_nvector_or_fail(
        lambda: N_VNewWithMemHelp_Cuda(5, use_managed_mem, mem_helper, sunctx)
    )
    _check_cuda_nvector_const(nvec, 4.0)


def _torch_dtype():
    import torch

    if sunrealtype == np.float32:
        return torch.float32
    if sunrealtype == np.float64:
        return torch.float64
    return torch.longdouble


def _jax_cuda_device_or_skip(jax):
    try:
        devices = [device for device in jax.devices() if device.platform in ("cuda", "gpu")]
    except RuntimeError as err:
        pytest.skip(f"JAX CUDA/GPU backend is not available: {err}")

    if not devices:
        pytest.skip("JAX CUDA/GPU backend is not available")
    return devices[0]


def _jax_dtype(jnp):
    if sunrealtype == np.float32:
        return jnp.float32
    if sunrealtype == np.float64:
        return jnp.float64
    return jnp.longdouble


@pytest.mark.skipif("N_VMake_Cuda" not in globals(), reason="CUDA bindings are not enabled")
def test_make_nvector_cuda_cupy_array(sunctx):
    cupy = pytest.importorskip("cupy")

    h_arr = np.zeros(5, dtype=sunrealtype)
    d_arr = cupy.arange(5, dtype=sunrealtype)

    nvec = N_VMake_Cuda(5, h_arr, d_arr, sunctx)
    h_arr_ref = weakref.ref(h_arr)
    d_arr_ref = weakref.ref(d_arr)
    del h_arr, d_arr
    gc.collect()
    assert h_arr_ref() is not None
    assert d_arr_ref() is not None

    N_VConst(3.0, nvec)

    assert_allclose(cupy.asnumpy(d_arr_ref()), 3.0)

    view = N_VGetCupyArray(nvec)
    assert_allclose(cupy.asnumpy(view), 3.0)

    with pytest.raises(TypeError):
        N_VGetCupyArray(nvec, device="cpu")

    del view, nvec
    gc.collect()
    assert h_arr_ref() is None
    assert d_arr_ref() is None


@pytest.mark.skipif("N_VMake_Cuda" not in globals(), reason="CUDA bindings are not enabled")
def test_make_nvector_cuda_torch_tensor(sunctx):
    torch = pytest.importorskip("torch")
    if not torch.cuda.is_available():
        pytest.skip("PyTorch CUDA is not available")

    h_arr = np.zeros(5, dtype=sunrealtype)
    d_arr = torch.arange(5, device="cuda", dtype=_torch_dtype())

    nvec = N_VMake_Cuda(5, h_arr, d_arr, sunctx)
    h_arr_ref = weakref.ref(h_arr)
    d_arr_ref = weakref.ref(d_arr)
    del h_arr, d_arr
    gc.collect()
    assert h_arr_ref() is not None
    assert d_arr_ref() is not None

    N_VConst(4.0, nvec)

    assert_allclose(d_arr_ref().cpu().numpy(), 4.0)

    view = N_VGetTorchTensor(nvec)
    assert_allclose(view.cpu().numpy(), 4.0)

    host_copy = view.cpu()
    assert_allclose(host_copy.numpy(), 4.0)
    assert_allclose(h_arr_ref(), 0.0)

    del view, host_copy, nvec
    gc.collect()
    assert h_arr_ref() is None
    assert d_arr_ref() is None


@pytest.mark.skipif("N_VMake_Cuda" not in globals(), reason="CUDA bindings are not enabled")
def test_make_nvector_cuda_jax_array(sunctx):
    jax = pytest.importorskip("jax")
    if np.dtype(sunrealtype) == np.dtype(np.float64):
        jax.config.update("jax_enable_x64", True)

    import jax.numpy as jnp

    device = _jax_cuda_device_or_skip(jax)
    h_arr = np.zeros(5, dtype=sunrealtype)
    d_arr = jax.device_put(jnp.arange(5, dtype=_jax_dtype(jnp)), device)
    d_arr.block_until_ready()

    nvec = N_VMake_Cuda(5, h_arr, d_arr, sunctx)
    h_arr_ref = weakref.ref(h_arr)
    d_arr_ref = weakref.ref(d_arr)
    del h_arr, d_arr
    gc.collect()
    assert h_arr_ref() is not None
    assert d_arr_ref() is not None

    assert int(N_VGetDeviceArrayPointer(nvec)) == d_arr_ref().unsafe_buffer_pointer()

    view = N_VGetJaxArray(nvec)
    assert_allclose(np.asarray(view), np.asarray(d_arr_ref()))

    host_copy = np.asarray(view)
    assert_allclose(host_copy, np.asarray(d_arr_ref()))

    N_VConst(6.0, nvec)

    updated_host_copy = np.asarray(N_VGetJaxArray(nvec))
    assert_allclose(updated_host_copy, 6.0)

    del view, host_copy, updated_host_copy, nvec
    gc.collect()
    assert h_arr_ref() is None
    assert d_arr_ref() is None


@pytest.mark.skipif("N_VNew_Cuda" not in globals(), reason="CUDA bindings are not enabled")
def test_nvector_array_helpers_cuda_reject_numpy(sunctx):
    nvec = _cuda_nvector_or_fail(lambda: N_VNew_Cuda(5, sunctx))

    with pytest.raises(RuntimeError, match="not supported for CUDA N_Vectors"):
        N_VGetNumpyArray(nvec)
    with pytest.raises(TypeError):
        N_VGetNumpyArray(nvec, device="cuda")


@pytest.mark.skipif("N_VNew_Cuda" not in globals(), reason="CUDA bindings are not enabled")
def test_set_nvector_cuda_jax_array(sunctx):
    jax = pytest.importorskip("jax")
    if np.dtype(sunrealtype) == np.dtype(np.float64):
        jax.config.update("jax_enable_x64", True)

    import jax.numpy as jnp

    device = _jax_cuda_device_or_skip(jax)
    nvec = _cuda_nvector_or_fail(lambda: N_VNew_Cuda(5, sunctx))
    original_pointer = int(N_VGetDeviceArrayPointer(nvec))
    array = jax.device_put(jnp.arange(5, dtype=_jax_dtype(jnp)), device)
    array.block_until_ready()

    N_VSetJaxArray(array, nvec)
    assert int(N_VGetDeviceArrayPointer(nvec)) == original_pointer
    assert_allclose(np.asarray(N_VGetJaxArray(nvec)), np.arange(5))

    N_VConst(6.0, nvec)
    assert_allclose(np.asarray(array), np.arange(5))

    cpu_array = jax.device_put(jnp.full(5, 2.0, dtype=_jax_dtype(jnp)), jax.devices("cpu")[0])
    N_VSetJaxArray(cpu_array, nvec)
    assert int(N_VGetDeviceArrayPointer(nvec)) == original_pointer
    assert_allclose(np.asarray(N_VGetJaxArray(nvec)), 2.0)


@pytest.mark.skipif("N_VMake_Cuda" not in globals(), reason="CUDA bindings are not enabled")
def test_set_nvector_cuda_jax_ref(sunctx):
    jax = pytest.importorskip("jax")
    if np.dtype(sunrealtype) == np.dtype(np.float64):
        jax.config.update("jax_enable_x64", True)

    import jax.numpy as jnp

    device = _jax_cuda_device_or_skip(jax)
    nvec = _cuda_nvector_or_fail(lambda: N_VNew_Cuda(5, sunctx))
    d_arr = jax.device_put(jnp.arange(5, dtype=_jax_dtype(jnp)), device)
    d_arr.block_until_ready()
    ref = jax.ref.new_ref(d_arr)
    ref[...].block_until_ready()
    cpu_ref = jax.ref.new_ref(
        jax.device_put(jnp.zeros(5, dtype=_jax_dtype(jnp)), jax.devices("cpu")[0])
    )

    with pytest.raises(RuntimeError, match="requires a JAX Ref on a CUDA device"):
        N_VSetJaxArray(cpu_ref, nvec, copy=False)

    ref_weak = weakref.ref(ref)
    N_VSetJaxArray(ref, nvec)
    assert int(N_VGetDeviceArrayPointer(nvec)) == ref.unsafe_buffer_pointer()

    with pytest.raises(RuntimeError, match="copy=True requires a jax.Array"):
        N_VSetJaxArray(ref, nvec, copy=True)

    del d_arr, ref
    gc.collect()
    assert ref_weak() is not None

    N_VConst(7.0, nvec)

    N_VGetJaxArray(nvec).block_until_ready()
    assert_allclose(np.asarray(ref_weak()[...]), 7.0)
    host_copy = np.asarray(N_VGetJaxArray(nvec))
    assert_allclose(host_copy, 7.0)

    del host_copy, nvec
    gc.collect()
    assert ref_weak() is None


# Test an operation that involves vector arrays
@pytest.mark.parametrize("vector_type", ["serial"])
def test_nvlinearcombination(vector_type, sunctx):
    if vector_type == "serial":
        nvec1 = N_VNew_Serial(5, sunctx)
        nvec2 = N_VNew_Serial(5, sunctx)
    else:
        raise ValueError("Unknown vector type")

    arr1 = N_VGetArrayPointer(nvec1)
    arr1[:] = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=sunrealtype)

    arr2 = N_VGetArrayPointer(nvec2)
    arr2[:] = np.array([10.0, 20.0, 30.0, 40.0, 50.0], dtype=sunrealtype)

    c = np.array([1.0, 0.1], dtype=sunrealtype)
    X = [nvec1, nvec2]

    z = N_VNew_Serial(5, sunctx)
    N_VConst(0.0, z)

    N_VLinearCombination(2, c, X, z)

    assert_allclose(N_VGetArrayPointer(z), [2.0, 4.0, 6.0, 8.0, 10.0])


def test_nvscaleaddmultivectorarray_serial(sunctx):
    nvec = 2
    nsum = 2
    length = 3

    # c_1d shape (nsum,)
    c_1d = np.array([2.0, 3.0], dtype=sunrealtype)

    # X_1d shape (nvec,)
    X_1d = [N_VNew_Serial(length, sunctx) for _ in range(nvec)]

    for i, x in enumerate(X_1d):
        N_VConst(float(i + 1), x)

    # Y_2d shape (nsum, nvec)
    Y_2d = [[N_VNew_Serial(length, sunctx) for _ in range(nvec)] for _ in range(nsum)]
    for s in range(nsum):
        for v in range(nvec):
            N_VConst(float((s + 1) * 10 + v), Y_2d[s][v])

    # Z_2d shape (nsum, nvec)
    Z_2d = [[N_VNew_Serial(length, sunctx) for _ in range(nvec)] for _ in range(nsum)]

    err = N_VScaleAddMultiVectorArray(nvec, nsum, c_1d, X_1d, Y_2d, Z_2d)
    assert err == SUN_SUCCESS

    # Check Z_2d[s][v] = c_1d[s] * X_1d[v] + Y_2d[s][v]
    for s in range(nsum):
        for v in range(nvec):
            expected = c_1d[s] * N_VGetArrayPointer(X_1d[v]) + N_VGetArrayPointer(Y_2d[s][v])
            actual = N_VGetArrayPointer(Z_2d[s][v])
            assert_allclose(actual, expected)


def test_nvlinearcombinationvectorarray_serial(sunctx):
    nvec = 2
    nsum = 2
    length = 3

    # c_1d shape (nsum,)
    c_1d = np.array([2.0, 3.0], dtype=sunrealtype)

    # X_2d shape (nsum, nvec)
    X_2d = []
    for s in range(nsum):
        row = []
        for v in range(nvec):
            x = N_VNew_Serial(length, sunctx)
            N_VConst(float((s + 1) * 10 + v), x)
            row.append(x)
        X_2d.append(row)

    # Z_1d shape (nvec,)
    Z_1d = [N_VNew_Serial(length, sunctx) for _ in range(nvec)]

    err = N_VLinearCombinationVectorArray(nvec, nsum, c_1d, X_2d, Z_1d)
    assert err == SUN_SUCCESS

    # Check Z_1d[v] = sum_s c_1d[s] * X_2d[s][v]
    for v in range(nvec):
        expected = sum(c_1d[s] * N_VGetArrayPointer(X_2d[s][v]) for s in range(nsum))
        actual = N_VGetArrayPointer(Z_1d[v])
        assert_allclose(actual, expected)
