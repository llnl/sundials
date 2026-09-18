# CVODE examples

Examples are grouped by problem. Each problem directory contains one or more
variants named for their implementation language and optional TPLs. For
example, `cvAdvDiff/c-mpi-hypre` is the MPI/HYPRE C variant of the
advection-diffusion example. Each variant contains its source files, reference
output, input data, supporting scripts, and `CMakeLists.txt`.

```text
examples/cvode/
├── README.md
├── CMakeLists.txt
├── cvAdvDiff/
│   ├── c/
│   ├── c-lapack/
│   ├── c-mpi/
│   ├── c-mpi-hypre/
│   ├── c-mpi-petsc/
│   ├── c-openmp/
│   ├── c-openmp-device/
│   ├── cpp-mpi-superludist/
│   ├── cpp-raja/
│   ├── cuda/
│   ├── cuda-managed/
│   ├── fortran/
│   ├── hip/
│   └── sycl/
├── cvAnalytic/
│   ├── c/
│   └── fortran/
├── cvAnalyticSys/
│   ├── fortran/
│   └── fortran-klu/
├── cvBrusselator/
│   ├── cpp-ginkgo/
│   ├── cpp-kokkos/
│   ├── cpp-kokkos-2d/
│   ├── cpp-magma/
│   └── fortran/
├── cvDiag/
│   └── fortran-mpi/
├── cvDirectDemo/
│   └── c/
├── cvDisc/
│   └── c/
├── cvDiurnal/
│   ├── c/
│   ├── c-mpi/
│   ├── c-mpi-mpimanyvector/
│   └── fortran/
├── cvHeat2D/
│   ├── cpp/
│   ├── cpp-ginkgo/
│   ├── cpp-mpi/
│   └── cpp-mpi-hypre/
├── cvKPR/
│   ├── cpp/
│   └── cpp-ginkgo/
├── cvKrylovDemo/
│   └── c/
├── cvParticle/
│   └── c/
├── cvPendulum/
│   └── c/
├── cvPetsc/
│   └── c-mpi-petsc/
├── cvRoberts/
│   ├── c/
│   ├── c-klu/
│   ├── c-lapack/
│   ├── c-superlumt/
│   ├── cpp-onemkl/
│   ├── cuda-cusolversp/
│   ├── fortran/
│   ├── fortran-klu/
│   └── fortran-lapack/
├── cvRocket/
│   └── c/
└── cvVdp/
    └── c/
```

The complete variant tree is visible directly under each problem directory.
Common variant names include `c`, `c-mpi`, `c-mpi-hypre`, `c-klu`, `cpp`,
`cpp-ginkgo`, `cpp-mpi`, `cuda`, `hip`, `sycl`, and `fortran`. CMake feature
guards select variants according to the enabled SUNDIALS modules and
third-party libraries; a variant that is not available in a build is skipped.
