#!/usr/bin/env bash
#SBATCH --job-name=diffusion-2d-py-cpp
#SBATCH --partition=pbatch
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:2
#SBATCH --time=04:00:00

set -euo pipefail

# Submit from the repository root, or set SUNDIALS_REPO_DIR.
repo_dir=${SUNDIALS_REPO_DIR:-${SLURM_SUBMIT_DIR:?}}
script_dir="$repo_dir/scripts/matrix"
python_example="$repo_dir/benchmarks/diffusion_2D/python/diffusion_2D_torch_mpi.py"
cpp_example="$repo_dir/build_diffusion_2d_mpicuda/benchmarks/diffusion_2D/mpi_gpu/arkode_diffusion_2D_mpicuda"
python_bin="$repo_dir/.venv-mpi-torch/bin/python"
repeats=${REPEATS:-5}
tf=${TF:-1}

if [[ ! $repeats =~ ^[1-9][0-9]*$ ]]; then
  echo "REPEATS must be a positive integer" >&2
  exit 1
fi
if [[ ! -x $cpp_example || ! -x $python_bin || ! -f $python_example ]]; then
  echo "Build the C++ benchmark and install sundials4py before submitting" >&2
  exit 1
fi

# shellcheck source=scripts/matrix/mpi_cuda_torch.profile
source "$script_dir/mpi_cuda_torch.profile"
"$python_bin" -c 'import mpi4py.MPI, torch, triton; import sundials4py.arkode; import sundials4py.core as sun; assert hasattr(sun, "N_VMake_MPIPlusX") and hasattr(sun, "N_VMake_Cuda")'

results_dir=${RESUME_DIR:-"$repo_dir/diffusion_2d_compare_${SLURM_JOB_ID}"}
csv="$results_dir/timings.csv"
stats_csv="$results_dir/integrator_stats.csv"
if [[ -n ${RESUME_DIR:-} ]]; then
  if [[ ! -f $csv || ! -f $stats_csv ]]; then
    echo "RESUME_DIR must contain timings.csv and integrator_stats.csv" >&2
    exit 1
  fi
  # An interrupted write to just one CSV needs repair before resuming.
  "$python_bin" - "$csv" "$stats_csv" <<'PY'
import csv
import sys

def keys(path):
    with open(path, newline="") as file:
        rows = list(csv.DictReader(file))
    result = [(row["gpus"], row["repeat"], row["rhs_mode"]) for row in rows]
    if len(result) != len(set(result)):
        raise RuntimeError(f"Duplicate benchmark rows in {path}")
    return set(result)

if keys(sys.argv[1]) != keys(sys.argv[2]):
    raise RuntimeError("Timing and integrator statistics rows differ")
PY
else
  mkdir -p "$results_dir"
  echo 'gpus,nodes,nx,ny,rank0_unknowns,other_rank_unknowns,repeat,rhs_mode,cpp_seconds,python_seconds,cpp_final_error,python_final_error' > "$csv"
  echo 'gpus,repeat,rhs_mode,cpp_steps,python_steps,cpp_step_attempts,python_step_attempts,cpp_nls_iters,python_nls_iters,cpp_ls_iters,python_ls_iters,cpp_prec_setups,python_prec_setups,cpp_prec_solves,python_prec_solves' > "$stats_csv"
fi

completed() {
  awk -F, -v g="$1" -v r="$2" -v m="$3" \
    'NR > 1 && $1 == g && $7 == r && $8 == m {found = 1} END {exit !found}' "$csv"
}

# nx=321 and dy=dx=1/320. Grow the y domain by one unit per rank.
# Rank 0 owns 321*321=103041 values; every other rank owns 321*320=102720.
nx=321
rows=320
nout=20

for gpus in 1 2 4 8; do
  nodes=$(( (gpus + 1) / 2 ))
  step_gpus=$(( gpus == 1 ? 1 : 2 ))
  ny=$(( rows * gpus + 1 ))
  echo "${gpus} GPU(s), ${nodes} node(s), grid ${nx}x${ny}, local values 103041/102720"

  for (( repeat=1; repeat<=repeats; repeat++ )); do
    if completed "$gpus" "$repeat" inplace && completed "$gpus" "$repeat" triton; then
      echo "Reusing completed $gpus GPU repeat $repeat"
      continue
    fi
    cpp_log="$results_dir/${gpus}gpu_repeat${repeat}_cpp.log"

    if [[ ! -s $cpp_log ]] || ! grep -q '^SOLVE_SECONDS=' "$cpp_log" ||
       ! grep -q '^FINAL_ERROR=' "$cpp_log" || ! grep -q '^Prec solves[[:space:]]*=' "$cpp_log"; then
      srun --nodes="$nodes" --ntasks="$gpus" --ntasks-per-node=2 \
        --gres="gpu:$step_gpus" --gpu-bind=single:1 --kill-on-bad-exit=1 \
        "$cpp_example" --npx 1 --npy "$gpus" --nx "$nx" --ny "$ny" \
        --xu 1 --yu "$gpus" --kx 1 --ky 1 --tf "$tf" --nout "$nout" \
        --rtol 1e-5 --atol 1e-10 --order 3 --controller I --ls cg \
        --liniters 20 --msbp 0 --epslin 0 --output 1 --host-staging \
        | tee "$cpp_log"
    else
      echo "Reusing completed C++ log: $cpp_log"
    fi

    cpp_seconds=$(awk -F= '/^SOLVE_SECONDS=/{print $2}' "$cpp_log")
    cpp_error=$(awk -F= '/^FINAL_ERROR=/{print $2}' "$cpp_log")
    cpp_steps=$(awk '/^Steps[[:space:]]+=/{print $NF}' "$cpp_log")
    cpp_attempts=$(awk '/^Step attempts[[:space:]]+=/{print $NF}' "$cpp_log")
    cpp_nls_iters=$(awk '/^NLS iters[[:space:]]+=/{print $NF}' "$cpp_log")
    cpp_ls_iters=$(awk '/^LS iters[[:space:]]+=/{print $NF}' "$cpp_log")
    cpp_prec_setups=$(awk '/^Prec setup evals[[:space:]]+=/{print $NF}' "$cpp_log")
    cpp_prec_solves=$(awk '/^Prec solves[[:space:]]+=/{print $NF}' "$cpp_log")
    if [[ -z $cpp_seconds || -z $cpp_error || -z $cpp_steps || -z $cpp_attempts ||
          -z $cpp_nls_iters || -z $cpp_ls_iters || -z $cpp_prec_setups || -z $cpp_prec_solves ]]; then
      echo "Missing C++ timing, error, or statistics in repeat $repeat at $gpus GPUs" >&2
      exit 1
    fi
    awk -v error="$cpp_error" 'BEGIN { exit !(error+0 <= 0.01) }' || {
      echo "C++ solution error exceeds 0.01 at $gpus GPUs" >&2
      exit 1
    }

    # Alternate Python order between repeats to reduce ordering bias.
    rhs_modes=(inplace triton)
    if (( repeat % 2 == 0 )); then rhs_modes=(triton inplace); fi
    for rhs_mode in "${rhs_modes[@]}"; do
      if completed "$gpus" "$repeat" "$rhs_mode"; then
        echo "Reusing completed $gpus GPU repeat $repeat $rhs_mode"
        continue
      fi
      python_log="$results_dir/${gpus}gpu_repeat${repeat}_python_${rhs_mode}.log"
      if [[ ! -s $python_log ]] || ! grep -q '^SOLVE_SECONDS=' "$python_log" ||
         ! grep -q '^FINAL_ERROR=' "$python_log" || ! grep -q '^PREC_SOLVES=' "$python_log"; then
        # Keep Triton's compiled kernels on node-local storage, isolated by rank.
        # The previous 4-GPU run hit a stale file handle in Triton's cache.
        # shellcheck disable=SC2016
        srun --nodes="$nodes" --ntasks="$gpus" --ntasks-per-node=2 \
          --gres="gpu:$step_gpus" --gpu-bind=single:1 --kill-on-bad-exit=1 \
          bash -c 'export TRITON_CACHE_DIR="${SLURM_TMPDIR:-/tmp}/triton-${SLURM_JOB_ID}-${SLURM_PROCID}"; mkdir -p "$TRITON_CACHE_DIR"; exec "$@"' bash \
          "$python_bin" "$python_example" --device cuda --nx "$nx" \
          --rows-per-rank "$rows" --kx 1 --ky 1 --tf "$tf" --nout "$nout" \
          --rhs-mode "$rhs_mode" | tee "$python_log"
      else
        echo "Reusing completed Python log: $python_log"
      fi

      python_seconds=$(awk -F= '/^SOLVE_SECONDS=/{print $2}' "$python_log")
      python_error=$(awk -F= '/^FINAL_ERROR=/{print $2}' "$python_log")
      python_steps=$(awk -F= '$1 == "STEPS" {print $2}' "$python_log")
      python_attempts=$(awk -F= '$1 == "STEP_ATTEMPTS" {print $2}' "$python_log")
      python_nls_iters=$(awk -F= '$1 == "NLS_ITERS" {print $2}' "$python_log")
      python_ls_iters=$(awk -F= '$1 == "LS_ITERS" {print $2}' "$python_log")
      python_prec_setups=$(awk -F= '$1 == "PREC_SETUPS" {print $2}' "$python_log")
      python_prec_solves=$(awk -F= '$1 == "PREC_SOLVES" {print $2}' "$python_log")
      if [[ -z $python_seconds || -z $python_error || -z $python_steps || -z $python_attempts ||
            -z $python_nls_iters || -z $python_ls_iters || -z $python_prec_setups || -z $python_prec_solves ]]; then
        echo "Missing Python $rhs_mode timing, error, or statistics in repeat $repeat at $gpus GPUs" >&2
        exit 1
      fi
      awk -v error="$python_error" 'BEGIN { exit !(error+0 <= 0.01) }' || {
        echo "Python $rhs_mode solution error exceeds 0.01 at $gpus GPUs" >&2
        exit 1
      }
      echo "$gpus,$nodes,$nx,$ny,103041,102720,$repeat,$rhs_mode,$cpp_seconds,$python_seconds,$cpp_error,$python_error" >> "$csv"
      echo "$gpus,$repeat,$rhs_mode,$cpp_steps,$python_steps,$cpp_attempts,$python_attempts,$cpp_nls_iters,$python_nls_iters,$cpp_ls_iters,$python_ls_iters,$cpp_prec_setups,$python_prec_setups,$cpp_prec_solves,$python_prec_solves" >> "$stats_csv"
    done
  done
done

summary="$results_dir/summary.csv"
"$python_bin" - "$csv" "$summary" "$repeats" <<'PY'
import csv
import statistics
import sys
from collections import defaultdict

timings_path, summary_path, expected_repeats = sys.argv[1], sys.argv[2], int(sys.argv[3])
groups = defaultdict(lambda: {"cpp": [], "python": []})
with open(timings_path, newline="") as timing_file:
    for row in csv.DictReader(timing_file):
        group = groups[(int(row["gpus"]), row["rhs_mode"])]
        group["cpp"].append(float(row["cpp_seconds"]))
        group["python"].append(float(row["python_seconds"]))

with open(summary_path, "w", newline="") as summary_file:
    writer = csv.writer(summary_file)
    writer.writerow(
        (
            "gpus", "rhs_mode", "repeats", "cpp_median_seconds",
            "python_median_seconds", "cpp_mean_seconds", "python_mean_seconds",
            "python_over_cpp_median", "cpp_weak_efficiency_median",
            "python_weak_efficiency_median",
        )
    )
    for gpus in (1, 2, 4, 8):
        for mode in ("inplace", "triton"):
            group = groups[(gpus, mode)]
            baseline = groups[(1, mode)]
            if len(group["cpp"]) != expected_repeats:
                raise RuntimeError(f"Expected {expected_repeats} runs for {gpus} GPUs, {mode}")
            cpp_median = statistics.median(group["cpp"])
            python_median = statistics.median(group["python"])
            writer.writerow(
                (
                    gpus, mode, expected_repeats, f"{cpp_median:.6f}",
                    f"{python_median:.6f}", f"{statistics.mean(group['cpp']):.6f}",
                    f"{statistics.mean(group['python']):.6f}",
                    f"{python_median / cpp_median:.4f}",
                    f"{statistics.median(baseline['cpp']) / cpp_median:.4f}",
                    f"{statistics.median(baseline['python']) / python_median:.4f}",
                )
            )
PY

echo "Paired timings: $csv"
echo "Integrator statistics: $stats_csv"
echo "Summary: $summary"
cat "$summary"
