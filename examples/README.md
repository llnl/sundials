# SUNDIALS examples

Examples are grouped first by SUNDIALS package and then by problem. Each
problem directory contains one or more variants named for their implementation
language and optional third-party libraries. For example,
`cvode/cv_adv_diff/c-mpi-hypre` is the MPI/HYPRE C variant of the CVODE
advection-diffusion example.

The directory layout is:

```text
examples/
└── <package>/
    └── <problem>/
        └── <language-and-TPL-variant>/
```

Browse the package and problem directories to find an example, then choose a
variant based on its implementation language and optional dependencies. Each
variant contains its source files, reference output, input data, supporting
scripts, and a `CMakeLists.txt` for building it against an installed SUNDIALS.
