#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd -- "$script_dir/../.." && pwd)
# shellcheck source=scripts/matrix/arkode_torch_cuda_env.sh
source "$script_dir/arkode_torch_cuda_env.sh"

venv_dir="$repo_dir/build_matrix_cuda/venv"
python3 -m venv "$venv_dir"
venv_python="$venv_dir/bin/python"

"$venv_python" -m pip install --upgrade pip
"$venv_python" -m pip install torch --index-url https://download.pytorch.org/whl/cu126

# The Matrix NVIDIA nodes have H100 GPUs (compute capability 9.0). Set the
# architecture explicitly because the login nodes do not have visible GPUs.
export CMAKE_ARGS="-DSUNDIALS_ENABLE_CUDA=ON -DCMAKE_CUDA_ARCHITECTURES=90 ${CMAKE_ARGS:-}"
"$venv_python" -m pip install -e "$repo_dir"

"$venv_python" - <<'PY'
import sundials4py.core as sun
import torch

assert hasattr(sun, "N_VMake_Cuda"), "sundials4py was built without CUDA"
assert torch.version.cuda is not None, "installed PyTorch has no CUDA support"
print(f"Installed CUDA sundials4py and PyTorch {torch.__version__}")
PY

echo "Submit from $repo_dir with: sbatch scripts/matrix/run_arkode_torch_cuda.sbatch"
