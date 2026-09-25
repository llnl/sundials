# CVODE examples

Examples are grouped by problem. Each problem directory contains one or more
variants named for their implementation language and optional TPLs. For
example, `cvAdvDiff/c-mpi-hypre` is the MPI/HYPRE C variant of the
advection-diffusion example. Each variant contains its source files, reference
output, input data, supporting scripts, and `CMakeLists.txt`.

The directory layout is:

```text
examples/cvode/
├── CMakeLists.txt
└── <problem>/
    └── <language-and-TPL-variant>/
```

Browse the problem directories to find an example, then choose a variant
based on its implementation language and optional TPLs. Common variant names
include `c`, `c-mpi`, `c-mpi-hypre`, `c-klu`, `cpp`, `cpp-ginkgo`, `cpp-mpi`,
`cuda`, `hip`, `sycl`, and `fortran`. Each variant contains its source files,
reference output, input data, supporting scripts, and `CMakeLists.txt`.
CMake feature guards select variants according to the enabled SUNDIALS
modules and third-party libraries; a variant that is not available in a build
is skipped.
