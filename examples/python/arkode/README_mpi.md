# Run MPI heat examples

Build sundials4py with MPI support. For the JAX and Torch variants, also enable
CUDA. Install mpi4py against the same MPI library that SUNDIALS uses.

```sh
python -m pip install mpi4py 'scikit-build-core>=1.1.2' 'nanobind>=2.9.2'
python -m pip install --no-build-isolation -e . \
	-Ccmake.define.SUNDIALS_ENABLE_MPI=ON \
	-Ccmake.define.SUNDIALS_ENABLE_CUDA=ON \
	-Ccmake.define.SUNDIALS_ENABLE_CUDA_EXAMPLES=OFF
python -m pip install 'jax[cuda13]' torch
```

Select the JAX CUDA extra that matches your CUDA installation. Use a CUDA-enabled
Torch distribution. If MPI libraries are outside the system search path, add
their directory to `LD_LIBRARY_PATH` before launching Python.

From the repository root, run the CPU example:

```sh
mpiexec -n 2 python -m mpi4py examples/python/arkode/ark_heat1D_mpi.py
```

Run the same problem with CUDA local vectors and JAX or Torch callbacks:

```sh
export XLA_PYTHON_CLIENT_PREALLOCATE=false
mpiexec -n 2 python -m mpi4py examples/python/arkode/ark_heat1D_mpi_jax.py
mpiexec -n 2 python -m mpi4py examples/python/arkode/ark_heat1D_mpi_torch.py
```

Each example prints solver work, analytic error, warmed integration time, setup
time, and warmup time. All ranks use visible GPU 0. On multi-GPU systems, use your
launcher to set `CUDA_VISIBLE_DEVICES` for each rank. On a single GPU, ranks share
that GPU. JAX preallocation must be disabled before import when ranks share it.

For ROCm, build sundials4py from source with `SUNDIALS_ENABLE_MPI=ON` and
`SUNDIALS_ENABLE_HIP=ON`, and install mpi4py against the same MPI library.
The separate `sundials4py-rocm` wheel is built without MPI and takes precedence
over the source extension when installed, so remove it from this environment.
Use the following source build, with `ROCM_PATH` and `HIP_PATH` set if CMake
cannot find the ROCm installation:

```sh
python -m pip uninstall sundials4py-rocm
python -m pip install --no-build-isolation -e . \
	-Ccmake.define.SUNDIALS_ENABLE_MPI=ON \
	-Ccmake.define.SUNDIALS_ENABLE_HIP=ON
```

Install ROCm builds of PyTorch and JAX, then run:

```sh
export XLA_PYTHON_CLIENT_PREALLOCATE=false
mpiexec -n 2 python -m mpi4py examples/python/arkode/ark_heat1D_mpi_jax_hip.py
mpiexec -n 2 python -m mpi4py examples/python/arkode/ark_heat1D_mpi_torch_hip.py
```

ROCm PyTorch uses the `torch.cuda` API to select and synchronize HIP devices.
Set `ROCR_VISIBLE_DEVICES` per rank when assigning separate GPUs. The HIP
examples also stage halo rows through host memory, so HIP-aware MPI is not
required.

To increase vector work without increasing spatial stiffness, pass `--fields`.
Fields are independent heat equations, not a coupled PDE. Pass `--multimode` to
initialize two sine modes with field-dependent amplitudes. The default remains a
single sine mode on 201 interior grid rows.

For a stabilized explicit diffusion solve, pass `--solver sts`:

```sh
mpiexec -n 2 python -m mpi4py examples/python/arkode/ark_heat1D_mpi_jax.py \
	--fields 256 --multimode --solver sts
```

The default implicit ARKStep solver uses SPGMR, classical Gram-Schmidt, and an
analytic Jacobian-vector callback. It has no preconditioner. The STS option uses
LSRKStep RKC2 with an analytic spectral bound and avoids Krylov solves.

JAX reads local CUDA views and copies immutable callback results back through
`N_VSetJaxArray`. Torch writes mutable output views and synchronizes before
returning to SUNDIALS. Torch compiles the stencil by default. Pass `--eager` to
disable compilation. Both versions stage only halo rows through host buffers;
they do not require CUDA-aware MPI. Keep the GPU scripts beside their CPU sibling,
which owns the shared problem, MPI exchange, solver configuration, and validation.

## Compare performance

Run the rotating benchmark from the repository root:

```sh
python scripts/benchmark_sundials4py_mpi.py \
	--ranks 1 2 --fields 1 256 4096 --repeats 5 \
	--output /tmp/mpi-heat-results.jsonl
```

The script launches backends in rotating order with one complete warmup and one
measured solve per process. It uses STS and the multimode initial condition for
all backends. It sets CPU math-library threads to one and disables JAX
preallocation. The launcher arguments use Open MPI syntax.

`solve_seconds` measures `ARKodeEvolve` through completion of accelerator work,
using the maximum elapsed time across ranks. It includes callbacks, halo transfers,
MPI reductions, and solver vector operations. Vector allocation, initialization,
solver setup, imports, JIT compilation during warmup, verification, and printing
are outside that interval. `setup_seconds` reports fresh vector and solver setup.
`process_seconds` includes process startup, imports, preparation, warmup, solve,
verification, printing, and teardown as observed by the benchmark driver.

The analytic check uses the exact solution of the spatial discretization, with
an absolute maximum-error limit of `2e-5`. The JAX timing boundary explicitly waits
on the final SUNDIALS vector's CUDA stream, not only JAX's effect queue.

The JSONL file retains every sample, numerical settings, backend versions,
correctness checks, and solver counts. The script rejects failed, incomplete,
non-double-precision, or incorrect solves. Inspect work counts before comparing
timings. A single shared GPU does not measure multi-GPU scaling. Desktop GPU use
and other machine activity can affect these measurements.

See [MPI heat performance on a shared GPU](mpi_heat_performance.md) for one
workstation's measured results and profiling limitations.
