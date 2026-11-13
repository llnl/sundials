#!/bin/bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd -- "$script_dir/../.." && pwd)
# shellcheck source=scripts/matrix/mpi_cuda_torch.profile
source "$script_dir/mpi_cuda_torch.profile"

# Keep this MPI build separate from the existing CUDA-only editable build.
# Suppress deprecation warnings from C, C++, and CUDA host code in this build.
export CMAKE_ARGS="-DSUNDIALS_ENABLE_CUDA=ON -DSUNDIALS_ENABLE_MPI=ON -DCMAKE_CUDA_ARCHITECTURES=90 -DCMAKE_C_FLAGS=-Wno-deprecated-declarations -DCMAKE_CXX_FLAGS=-Wno-deprecated-declarations -DCMAKE_CUDA_FLAGS=-Xcompiler=-Wno-deprecated-declarations ${CMAKE_ARGS:-}"
python -m pip install -e "$repo_dir" \
    --config-settings="build-dir=$repo_dir/build_sundials4py_mpi_cuda_torch/{wheel_tag}"

python -c 'import mpi4py.MPI; import sundials4py.arkode; import sundials4py.core as sun; import torch; assert hasattr(sun, "N_VMake_Cuda"); assert torch.version.cuda is not None; print("sundials4py MPI/CUDA Torch environment installed")'
