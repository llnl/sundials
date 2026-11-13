#!/bin/bash
set -euo pipefail

# Run inside an allocation with at least two NVIDIA GPUs.
srun --ntasks=2 --gpu-bind=single:1 --kill-on-bad-exit=1 python -c '
from mpi4py import MPI
import torch

comm = MPI.COMM_WORLD
assert comm.size == 2, f"Expected two MPI ranks, got {comm.size}"
assert torch.version.cuda is not None, "PyTorch was built without CUDA"
assert torch.cuda.is_available(), f"Rank {comm.rank} cannot access CUDA"

local_comm = comm.Split_type(MPI.COMM_TYPE_SHARED)
local_rank = local_comm.rank
local_comm.Free()
visible = torch.cuda.device_count()
device_id = 0 if visible == 1 else local_rank
assert device_id < visible, f"Rank {comm.rank} has no assigned GPU"
torch.cuda.set_device(device_id)
value = torch.tensor([comm.rank + 1], device="cuda")
assert comm.allreduce(value.item(), op=MPI.SUM) == 3

print(f"Rank {comm.rank}: CUDA device {torch.cuda.current_device()}, tensor {value.item()}")
if comm.rank == 0:
    print(f"MPI: {MPI.Get_library_version().splitlines()[0]}")
    print(f"PyTorch: {torch.__version__}, built for CUDA {torch.version.cuda}")
'
