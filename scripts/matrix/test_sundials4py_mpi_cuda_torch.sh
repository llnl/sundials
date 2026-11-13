#!/bin/bash
set -euo pipefail

# Run in a two-GPU allocation with the MPI/CUDA Torch profile environment.
script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=scripts/matrix/mpi_cuda_torch.profile
source "$script_dir/mpi_cuda_torch.profile"
srun --ntasks=2 --gpu-bind=single:1 --kill-on-bad-exit=1 python -c '
from mpi4py import MPI
import numpy as np
import torch
import sundials4py.arkode
import sundials4py.core as sun

comm = MPI.COMM_WORLD
assert comm.size == 2, f"Expected two MPI ranks, got {comm.size}"
assert hasattr(sun, "N_VMake_Cuda"), "sundials4py has no CUDA N_Vector"
assert torch.cuda.is_available(), f"Rank {comm.rank} cannot access CUDA"
local = comm.Split_type(MPI.COMM_TYPE_SHARED)
device = 0 if torch.cuda.device_count() == 1 else local.rank
local.Free()
torch.cuda.set_device(device)

status, ctx = sun.SUNContext_Create(sun.SUN_COMM_NULL)
assert status == sun.SUN_SUCCESS
dtype = torch.float32 if np.dtype(sun.sunrealtype) == np.dtype(np.float32) else torch.float64
host = np.zeros(3, dtype=sun.sunrealtype)
data = torch.zeros(3, dtype=dtype, device="cuda")
vec = sun.N_VMake_Cuda(3, host, data, ctx)
sun.N_VConst(float(comm.rank + 1), vec)
assert torch.all(sun.N_VGetTorchTensor(vec) == comm.rank + 1).item()
assert comm.allreduce(data.sum().item(), op=MPI.SUM) == 9
if comm.rank == 0:
    print("sundials4py + mpi4py + CUDA Torch: OK on 2 ranks")
'
