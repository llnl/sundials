#!/usr/bin/env bash
#SBATCH --job-name=arkode-torch-cuda
#SBATCH --partition=pdebug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --time=00:30:00

set -euo pipefail

# Slurm runs a copy of this script from its spool directory, so use the
# directory from which sbatch was called rather than BASH_SOURCE.
repo_dir=${SUNDIALS_REPO_DIR:-${SLURM_SUBMIT_DIR:?Submit this job with sbatch}}
script_dir="$repo_dir/scripts/matrix"
if [[ ! -f "$repo_dir/examples/python/arkode/ark_heat1D_torch.py" ]]; then
  echo "SUNDIALS_REPO_DIR (or the sbatch directory) must be the repository root" >&2
  exit 1
fi
# shellcheck source=scripts/matrix/arkode_torch_cuda_env.sh
source "$script_dir/arkode_torch_cuda_env.sh"
venv_python="$repo_dir/build_matrix_cuda/venv/bin/python"
if [[ ! -x "$venv_python" ]]; then
  echo "Run scripts/matrix/setup_arkode_torch_cuda.sh before submitting" >&2
  exit 1
fi

srun "$venv_python" "$repo_dir/examples/python/arkode/ark_heat1D_torch.py" --device cuda "$@"
