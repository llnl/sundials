#!/usr/bin/env bash
# Matrix modules for the MPI + PyTorch CUDA ARKODE example.

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=scripts/matrix/arkode_torch_cuda_env.sh
source "$script_dir/arkode_torch_cuda_env.sh"
module load mvapich2/2.3.7
