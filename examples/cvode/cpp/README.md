# CVODE C++ examples

Each subdirectory contains one CVODE C++ example, its answer files, and any
supporting files needed by that example.

The examples include:

- `cv_heat2D` and `cv_kpr`: serial CVODE examples.
- `cv_heat2D_p`: MPI parallel heat equation example.
- `cv_heat2D_hypre_ls` and `cv_heat2D_hypre_pfmg`: MPI/HYPRE examples.
- `cvAdvDiff_sludist`: MPI/SuperLU_DIST example.
- `cv_bruss_batched_magma`: MAGMA batched solver example.
- `cv_bruss_batched_ginkgo`, `cv_heat2D_ginkgo`, and `cv_kpr_ginkgo`:
  Ginkgo backend examples.
- `cv_bruss_batched_kokkos` and `cv_bruss_batched_kokkos_2D`: Kokkos backend
  examples.
- `cvRoberts_blockdiag_onemkl`: oneMKL example.
- `cvAdvDiff_kry_raja`: RAJA example.

Examples requiring optional packages are built only when the corresponding
SUNDIALS feature is enabled.
