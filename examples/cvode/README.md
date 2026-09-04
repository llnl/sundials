# CVODE examples

CVODE examples are grouped by implementation language or programming model.
Each example has its own directory containing its source files, reference
output, input data, and example-specific scripts.

- `c/` contains C examples, including serial, MPI, OpenMP, hypre, and PETSc
  examples.
- `cpp/` contains portable C++ examples, including examples using Ginkgo,
  Kokkos, MAGMA, oneMKL, RAJA, and SuperLU_DIST.
- `cuda/` contains native CUDA examples.
- `hip/` contains native HIP examples.
- `sycl/` contains native SYCL examples.
- `fortran/` contains serial and MPI Fortran 2003 examples.

Execution models and optional dependencies are intentionally not additional
directory levels because an example can use more than one of them. CMake only
adds an example when its required SUNDIALS modules and third-party libraries
are enabled.
