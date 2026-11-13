#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
exec bash "$script_dir/setup_sundials4py_mpi_cuda_torch.sh"
