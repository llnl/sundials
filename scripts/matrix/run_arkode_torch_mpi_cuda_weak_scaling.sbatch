#!/usr/bin/env bash
#SBATCH --job-name=arkode-py-weak-scaling
#SBATCH --partition=pbatch
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:2
#SBATCH --time=02:00:00

set -euo pipefail

# Use two ranks and GPUs per node so every scale has the same rank placement.
repo_dir=${SUNDIALS_REPO_DIR:-${SLURM_SUBMIT_DIR:?Submit this job with sbatch}}
script_dir="$repo_dir/scripts/matrix"
example="$repo_dir/examples/python/arkode/ark_heat1D_torch_mpi.py"
venv_python="$repo_dir/.venv-mpi-torch/bin/python"
points_per_gpu=${POINTS_PER_GPU:-100000}
repeats=${REPEATS:-3}

if [[ ! -f "$example" || ! -x "$venv_python" ]]; then
  echo "Run scripts/matrix/setup_sundials4py_mpi_cuda_torch.sh before submitting" >&2
  exit 1
fi
if [[ ! $points_per_gpu =~ ^[1-9][0-9]*$ || ! $repeats =~ ^[1-9][0-9]*$ ]]; then
  echo "POINTS_PER_GPU and REPEATS must be positive integers" >&2
  exit 1
fi
if (( points_per_gpu < 100000 )); then
  echo "POINTS_PER_GPU must be at least 100000" >&2
  exit 1
fi

# shellcheck disable=SC1090,SC1091
source "$script_dir/mpi_cuda_torch.profile"

results_dir="${SLURM_SUBMIT_DIR}/weak_scaling_${SLURM_JOB_ID}"
mkdir -p "$results_dir"
csv="$results_dir/timings.csv"
echo 'gpus,nodes,fields,unknowns_per_gpu,total_unknowns,repeat,solve_seconds' > "$csv"

for gpus in 1 2 4; do
  nodes=$(( (gpus + 1) / 2 ))
  fields=$(( (points_per_gpu + 98) / 99 ))
  unknowns_per_gpu=$(( fields * 99 ))
  total_unknowns=$(( gpus * unknowns_per_gpu ))
  step_gpus=2
  if (( gpus == 1 )); then step_gpus=1; fi

  for (( repeat=1; repeat<=repeats; repeat++ )); do
    log="$results_dir/${gpus}gpus_repeat${repeat}.out"
    echo "Running $gpus GPUs on $nodes nodes: $fields fields and $unknowns_per_gpu unknowns/GPU, repeat $repeat/$repeats"
    srun --nodes="$nodes" --ntasks="$gpus" --ntasks-per-node=2 \
      --gres="gpu:$step_gpus" --gpu-bind=single:1 \
      "$venv_python" "$example" --device cuda --n 101 --fields "$fields" \
      --timing | tee "$log"
    seconds=$(awk -F= '/^SOLVE_SECONDS=/{print $2}' "$log")
    if [[ -z $seconds ]]; then
      echo "No solver timing found in $log" >&2
      exit 1
    fi
    echo "$gpus,$nodes,$fields,$unknowns_per_gpu,$total_unknowns,$repeat,$seconds" >> "$csv"
  done
done

summary="$results_dir/summary.csv"
awk -F, '
  NR > 1 { sum[$1] += $7; count[$1]++ }
  END {
    print "gpus,mean_solve_seconds,weak_scaling_efficiency"
    baseline = sum[1] / count[1]
    split("1 2 4", scales, " ")
    for (i = 1; i <= 3; i++) {
      gpus = scales[i]
      mean = sum[gpus] / count[gpus]
      printf "%d,%.6f,%.6f\n", gpus, mean, baseline / mean
    }
  }
' "$csv" > "$summary"

echo "Weak scaling timings: $csv"
echo "Weak scaling summary: $summary"
