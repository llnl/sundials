# -----------------------------------------------------------------
# SUNDIALS Copyright Start
# Copyright (c) 2025-2026, Lawrence Livermore National Security,
# University of Maryland Baltimore County, and the SUNDIALS contributors.
# All rights reserved.
#
# SPDX-License-Identifier: BSD-3-Clause
# SUNDIALS Copyright End
# -----------------------------------------------------------------

import struct

import triton
import triton.language as tl


# Triton treats Python float kernel arguments as float32. Pass the changing
# forcing coefficients as tagged integer words to retain sunrealtype precision.
def coefficient_words(value, double):
    bits = struct.unpack(
        "<Q" if double else "<I", struct.pack("<d" if double else "<f", value)
    )[0]
    tag = 1 << 33  # Keeps the launch argument type stable even when a word is zero.
    return (bits & 0xFFFFFFFF) | tag, ((bits >> 32) & 0xFFFFFFFF) | tag


@triton.jit
def _coefficient(low, high, DOUBLE: tl.constexpr):
    if DOUBLE:
        bits = (tl.cast(high, tl.uint32).to(tl.uint64) << 32) | tl.cast(
            low, tl.uint32
        ).to(tl.uint64)
        return tl.cast(bits, tl.float64, bitcast=True)
    else:
        return tl.cast(tl.cast(low, tl.uint32), tl.float32, bitcast=True)


@triton.jit
def _interior(
    Y,
    YDOT,
    SHAPE,
    SPATIAL,
    A_LOW,
    A_HIGH,
    B_LOW,
    B_HIGH,
    NX: tl.constexpr,
    ROWS: tl.constexpr,
    CX: tl.constexpr,
    CY: tl.constexpr,
    CC: tl.constexpr,
    DOUBLE: tl.constexpr,
    SINGLE_RANK: tl.constexpr,
    BLOCK: tl.constexpr,
):
    offset = tl.program_id(0) * BLOCK + tl.arange(0, BLOCK)
    if SINGLE_RANK:
        valid = offset < ROWS * NX
        idx = offset
        row = offset // NX
    else:
        valid = offset < (ROWS - 2) * NX
        idx = offset + NX
    col = offset % NX
    active = valid & (col > 0) & (col < NX - 1)
    if SINGLE_RANK:
        active = active & (row > 0) & (row < ROWS - 1)

    center = tl.load(Y + idx, mask=active, other=0)
    west = tl.load(Y + idx - 1, mask=active, other=0)
    east = tl.load(Y + idx + 1, mask=active, other=0)
    south = tl.load(Y + idx - NX, mask=active, other=0)
    north = tl.load(Y + idx + NX, mask=active, other=0)
    shape = tl.load(SHAPE + idx, mask=active, other=0)
    spatial = tl.load(SPATIAL + idx, mask=active, other=0)
    a = _coefficient(A_LOW, A_HIGH, DOUBLE)
    b = _coefficient(B_LOW, B_HIGH, DOUBLE)
    result = CC * center + CX * (west + east) + CY * (south + north)
    result = result + a * shape + b * spatial
    tl.store(YDOT + idx, tl.where(active, result, 0.0), mask=valid)


@triton.jit
def _boundary(
    Y,
    YDOT,
    SHAPE,
    SPATIAL,
    SOUTH_HALO,
    NORTH_HALO,
    A_LOW,
    A_HIGH,
    B_LOW,
    B_HIGH,
    NX: tl.constexpr,
    ROWS: tl.constexpr,
    CX: tl.constexpr,
    CY: tl.constexpr,
    CC: tl.constexpr,
    HAS_SOUTH: tl.constexpr,
    HAS_NORTH: tl.constexpr,
    DOUBLE: tl.constexpr,
    BLOCK: tl.constexpr,
):
    side = tl.program_id(0)
    col = tl.program_id(1) * BLOCK + tl.arange(0, BLOCK)
    valid = col < NX
    has_halo = tl.where(side == 0, HAS_SOUTH, HAS_NORTH)
    active = valid & (col > 0) & (col < NX - 1) & has_halo
    row = tl.where(side == 0, 0, ROWS - 1)
    idx = row * NX + col

    center = tl.load(Y + idx, mask=active, other=0)
    west = tl.load(Y + idx - 1, mask=active, other=0)
    east = tl.load(Y + idx + 1, mask=active, other=0)
    inward = tl.load(Y + idx + tl.where(side == 0, NX, -NX), mask=active, other=0)
    south = tl.load(SOUTH_HALO + col, mask=active & (side == 0), other=0)
    north = tl.load(NORTH_HALO + col, mask=active & (side == 1), other=0)
    halo = tl.where(side == 0, south, north)
    shape = tl.load(SHAPE + idx, mask=active, other=0)
    spatial = tl.load(SPATIAL + idx, mask=active, other=0)
    a = _coefficient(A_LOW, A_HIGH, DOUBLE)
    b = _coefficient(B_LOW, B_HIGH, DOUBLE)
    result = CC * center + CX * (west + east) + CY * (inward + halo)
    result = result + a * shape + b * spatial
    tl.store(YDOT + idx, tl.where(active, result, 0.0), mask=valid)


def launch_interior(problem, y, ydot, words):
    single_rank = problem.south_recv is None and problem.north_recv is None
    count = (problem.local_rows if single_rank else problem.local_rows - 2) * problem.nx
    if count <= 0:
        return
    block = 256
    _interior[(triton.cdiv(count, block),)](
        y,
        ydot,
        problem.shape,
        problem.spatial_forcing,
        *words,
        problem.nx,
        problem.local_rows,
        problem.cx,
        problem.cy,
        problem.cc,
        problem.triton_double,
        single_rank,
        block,
    )


def launch_boundary(problem, y, ydot, words):
    if problem.south_recv is None and problem.north_recv is None:
        return
    block = 256
    _boundary[(2, triton.cdiv(problem.nx, block))](
        y,
        ydot,
        problem.shape,
        problem.spatial_forcing,
        problem.south_recv_device,
        problem.north_recv_device,
        *words,
        problem.nx,
        problem.local_rows,
        problem.cx,
        problem.cy,
        problem.cc,
        problem.south_recv is not None,
        problem.north_recv is not None,
        problem.triton_double,
        block,
    )
