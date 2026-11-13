# ARKODE PyTorch CUDA example on Matrix

From the repository root on a Matrix login node, install the CUDA build once:

```bash
bash scripts/matrix/setup_arkode_torch_cuda.sh
```

The setup uses Matrix's GCC 13.3.1, CUDA 12.6, Python 3.12, and CMake 3.30
modules. It installs a CUDA PyTorch wheel and an editable CUDA enabled
`sundials4py` build in `build_matrix_cuda/venv`.

Submit the example to a NVIDIA GPU node:

```bash
sbatch scripts/matrix/run_arkode_torch_cuda.sbatch
```

Pass example arguments after the job script, for example:

```bash
sbatch scripts/matrix/run_arkode_torch_cuda.sbatch --n 1001
```

The default job uses `pdebug` for 30 minutes. Slurm writes its output to
`slurm-<jobid>.out` in the submission directory. For a longer run, override
the partition and time at submission:

```bash
sbatch -p pbatch -t 02:00:00 scripts/matrix/run_arkode_torch_cuda.sbatch
```

## MPI + CUDA

Build `sundials4py` with MPI and CUDA into `.venv-mpi-torch`, the virtual
environment created by `mpi_cuda_torch.profile`:

```bash
bash scripts/matrix/setup_sundials4py_mpi_cuda_torch.sh
```

This build suppresses `-Wdeprecated-declarations` warnings in C, C++, and CUDA
host code.

Submit the distributed Python example from the repository root:

```bash
sbatch scripts/matrix/run_arkode_torch_mpi_cuda.sbatch
sbatch scripts/matrix/run_arkode_torch_mpi_cuda.sbatch --n 201
```

The job requests two MPI ranks and two NVIDIA GPUs on one node. Each rank uses
one GPU and exchanges heat equation boundary values with its neighbors through
`mpi4py`. The example uses fixed step explicit ARKODE integration to keep MPI
communications synchronized. `--n` is the global number of grid points and
must be at least the number of ranks plus two. Larger grids need more steps
because the explicit method's stable step size decreases with grid spacing.

## MPI + CUDA weak scaling

After building the MPI/CUDA environment, submit the weak scaling job from the
repository root:

```bash
sbatch scripts/matrix/run_arkode_torch_mpi_cuda_weak_scaling.sbatch
```

The job runs 1, 2, and 4 ranks, with one GPU per rank and at least 100,000
unknowns per GPU. It allocates two nodes with two GPUs per node. Each GPU
solves 1,011 independent copies (`--fields`) of the same 101-point heat
problem: 99 unknowns per copy, or 100,089 unknowns per GPU. Its grid spacing,
diffusivity, source, final time, and number of time steps stay fixed at every
scale. Each scale runs three times. Solver wall time is the maximum across
ranks. Initialization and exact-solution validation are outside the timer. Logs,
`timings.csv`, and `summary.csv` go into
`weak_scaling_<jobid>/` in the submission directory.
The summary reports mean solve time and weak scaling efficiency relative to
the one-GPU mean time. `POINTS_PER_GPU` is a minimum and is rounded up to a
whole number of 99-unknown fields.

Set `POINTS_PER_GPU` or `REPEATS` at submission to change the run, for
example:

```bash
sbatch --export=ALL,POINTS_PER_GPU=200000,REPEATS=5 scripts/matrix/run_arkode_torch_mpi_cuda_weak_scaling.sbatch
```

## CUDA Torch and mpi4py smoke tests

From the repository root, request two GPUs on one node, then source the
matching profile and run its test script:

```bash
salloc --partition=pdebug --nodes=1 --ntasks=2 --gres=gpu:2 --time=00:30:00
source scripts/matrix/cuda_torch.profile
bash scripts/matrix/test_torch_cuda.sh

source scripts/matrix/mpi_cuda_torch.profile
bash scripts/matrix/test_mpi_torch_cuda.sh
```

The MPI test starts two Slurm tasks and needs two GPUs in the allocation. The
MPI profile loads Matrix's GCC 13.3.1, MVAPICH2 2.3.7, and Python 3.12 modules;
the virtual environment inherits Matrix's CUDA PyTorch and mpi4py packages.

In the two-GPU allocation, check that `sundials4py` can use a CUDA Torch
tensor on both MPI ranks before running the heat example:

```bash
bash scripts/matrix/test_sundials4py_mpi_cuda_torch.sh
```
