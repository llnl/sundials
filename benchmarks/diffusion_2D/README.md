# Benchmark: 2D Diffusion

This benchmark problem implements a 2D diffusion equation using either MPI,
MPI + CUDA, or MPI + HIP parallelism. Note a GPU-aware MPI implementation is
required.

## Problem description

This code simulates the anisotropic 2D heat equation,

$$\frac{\partial u}{\partial t} = \nabla \cdot (D \nabla u) + b(t, \mathbf{x})$$

where $D$ is a diagonal matrix with entries $k_x$ and $k_y$. The system is
evolved for $t \in [0, t_f]$ on the rectangular domain
$(x,y) \equiv \mathbf{x} \in [\mathbf{0}, \mathbf{x}_{\text{max}}]^2$, with the
initial condition

$$u(0,\mathbf{x}) = \sin^2(\pi x) \sin^2(\pi y),$$

and stationary boundary conditions

$$\frac{\partial u}{\partial t}(t,0,y) = \frac{\partial u}{\partial t}(t,x_{\text{max}},y) = \frac{\partial u}{\partial t}(t,x,0) = \frac{\partial u}{\partial t}(t,x,y_{\text{max}}) = 0.$$

The source term is given by

$$b(t,\mathbf{x}) = -2 \pi \sin^2(\pi x) \sin^2(\pi y) \sin(\pi t) \cos(\pi t) - k_x 2 \pi^2 (\cos^2(\pi x) - \sin^2(\pi x)) \sin^2(\pi y) \cos^2(\pi t) - k_y 2 \pi^2 (\cos^2(\pi y) - \sin^2(\pi y)) \sin^2(\pi x) \cos^2(\pi t).$$

Under this setup, the problem has the analytical solution

$$u(t,\mathbf{x}) = \sin^2(\pi x) \sin^2(\pi y) \cos^2(\pi t).$$

Spatial derivatives are computed using second-order centered differences on a
uniform spatial grid. The problem can be evolved in time with ARKODE, CVODE, or
IDA. With ARKODE, an adaptive step diagonally implicit Runge-Kutta (DIRK) method
is applied. When using CVODE or IDA, adaptive order and step BDF methods are
used.

By default, the nonlinear system(s) in each time step are solved using an
inexact Newton method paired with a matrix-free CG linear solver and a Jacobi
preconditioner. A matrix-free GMRES linear solver may be selected at run time.
If SUNDIALS is built with the SuperLU_DIST interface enabled a modified Newton
method with SuperLU_DIST as the direct linear solver may also be selected at run
time.

## Options

Several command line options are available to change the problem parameters
as well as the integrator and solver options. A summary of the options are
listed below.

| Option                               | Description                                                                              | Default |
|:-------------------------------------|:-----------------------------------------------------------------------------------------|:--------|
| `--help`                             | Print the command line options and description                                           | --      |
| Problem Configuration Options        |                                                                                          |         |
| `--npx <int>`                        | Number of MPI tasks in the x-direction (0 forces MPI to decide)                          | 0       |
| `--npy <int>`                        | Number of MPI tasks in the y-direction (0 forces MPI to decide)                          | 0       |
| `--nx <int>`                         | Number of mesh points in the x-direction                                                 | 32      |
| `--ny <int>`                         | Number of mesh points in the y-direction                                                 | 32      |
| `--xu <sunrealtype>`                 | The domain upper bound in the x-direction $x_{\text{max}}$                               | 1.0     |
| `--yu <sunrealtype>`                 | The domain upper bound in the y-direction $y_{\text{max}}$                               | 1.0     |
| `--kx <sunrealtype>`                 | Diffusion coefficient in the x-direction $k_x$                                           | 1.0     |
| `--ky <sunrealtype>`                 | Diffusion coefficient in the y-direction $k_y$                                           | 1.0     |
| `--tf <sunrealtype>`                 | The final time `tf`                                                                      | 1.0     |
| `--noforcing`                        | Disable the forcing term                                                                 | Enabled |
| Output Options                       |                                                                                          |         |
| `--output <int>`                     | Output level: `0` no output, `1` output progress and stats, `2` write solution to disk   | 1       |
| `--nout <int>`                       | Number of output times                                                                   | 20      |
| Common Integrator and Solver Options |                                                                                          |         |
| `--rtol <sunrealtype>`               | Relative tolerance                                                                       | 1e-5    |
| `--atol <sunrealtype>`               | Absolute tolerance                                                                       | 1e-10   |
| `--maxsteps <int>`                   | Max number of steps between outputs (0 uses the integrator default)                      | 0       |
| `--onstep <int>`                     | Number of steps to run using `ONE_STEP` mode for debugging (0 uses `NORMAL` mode)        | 0       |
| `--ls <cg,gmres,sludist>`            | Linear solver: CG, GMRES, or SuperLU_DIST                                                | cg      |
| `--liniters <int>`                   | Number of linear iterations                                                              | 20      |
| `--epslin <sunrealtype>`             | Linear solve tolerance factor (0 uses the integrator default)                            | 0       |
| `--msbp <int>`                       | The linear solver setup frequency (CVODE and ARKODE only, 0 uses the integrator default) | 0       |
| Additional ARKODE Options            |                                                                                          |         |
| `--order <int>`                      | Methods order                                                                            | 3       |
| `--controller <I,PI,PID,...>`        | Error controller option                                                                  | 0       |
| `--nonlinear`                        | Treat the problem as nonlinearly implicit                                                | Linear  |

## Building

To build the benchmark executables SUNDIALS should be configured with ARKODE,
CVODE, or IDA enabled, MPI support turned on, and benchmarks enabled. If
SUNDIALS is configured with SuperLU_DIST enabled this linear solver can be
selected at run time and may utilizie OpenMP, CUDA, or ROCM (HIP) for on-node
parallelism. If SUNDIALS is configured with CUDA or HIP support enabled
additional executables utilizing CUDA and HIP will be built. See the SUNDIALS
installation guide for more details on configuring, building, and installing.

## Running

Based on the configuration, executables for each integrator and backend option
are built and installed in `<SUNDIALS_BENCHMARKS_INSTALL_PATH>/diffusion_2D`. The
executables follow the naming convention `<package>_diffusion_2D_<parallelism>`
where `<package>` is `arkode`, `cvode`, or `ida` and `<parallelism>` is `mpi` for
MPI only parallelism, `mpicuda` for MPI + CUDA, and `mpihip` for MPI + HIP.

**Note:** When using the SuperLU_DIST linear solver computations will be
offloaded to the GPU in the MPI only executables if CUDA or ROCM support is
enabled in SuperLU_DIST.

On Summit, with the default environment
```
  Compiler: xl/16.1.1-5
  MPI: spectrum-mpi/10.3.1.2-20200121
  CUDA: cuda/10.1.243
```
an example `jsrun` command using CUDA-aware MPI is
```
jsrun --smpiargs="-gpu" -n 2 -a 1 -c 1 -g 1 ./cvode_diffusion_2D_mpicuda
```

On Lassen, with the environment
```
  Compiler: gcc/8.3.1
  MPI: mvapich2/2021.05.28-cuda-11.1.1
  CUDA: cuda/11.1.1
```
an example `jsrun` command using CUDA-aware MPI is
```
jsrun -n 2 -a 1 -c 1 -g 1 ./cvode_diffusion_2D_mpicuda
```

On Crusher, with the environment
```
  Compiler: clang/14.0.2
  MPI: cray-mpich/8.1.17
  ROCM: rocm/5.2.0
```
an example `srun` command is
```
srun -N1 -n8 -c1 --gpus-per-node=8 --gpu-bind=closest ./cvode_diffusion_2D_mpi
```
