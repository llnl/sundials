#!/usr/bin/env bash
#SBATCH --job-name=arkode-torch-mpi-cuda
#SBATCH --partition=pdebug
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:2
#SBATCH --time=00:30:00

set -euo pipefail

# Slurm runs a spool copy, so BASH_SOURCE does not identify the repository.
repo_dir=${SUNDIALS_REPO_DIR:-${SLURM_SUBMIT_DIR:?Submit this job with sbatch}}
script_dir="$repo_dir/scripts/matrix"
example="$repo_dir/examples/python/arkode/ark_heat1D_torch_mpi.py"
if [[ ! -f "$example" ]]; then
  echo "SUNDIALS_REPO_DIR (or the sbatch directory) must be the repository root" >&2
  exit 1
fi
# shellcheck source=scripts/matrix/mpi_cuda_torch.profile
source "$script_dir/mpi_cuda_torch.profile"
venv_python="$repo_dir/.venv-mpi-torch/bin/python"
if [[ ! -x "$venv_python" ]] || ! "$venv_python" -c 'import mpi4py.MPI; import sundials4py.arkode; import sundials4py.core as sun; import torch; assert hasattr(sun, "N_VMake_Cuda")'; then
  echo "Run scripts/matrix/setup_sundials4py_mpi_cuda_torch.sh before submitting" >&2
  exit 1
fi

srun --gpu-bind=single:1 "$venv_python" "$example" --device cuda "$@"
