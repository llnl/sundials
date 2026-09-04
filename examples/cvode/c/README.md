# CVODE C examples

Each subdirectory contains one CVODE C example, its answer files, and any
supporting files needed by that example.

The examples include:

- `cvAdvDiff_bnd`, `cvAnalytic_mels`, `cvDirectDemo_ls`, `cvDisc_dns`,
  `cvDiurnal_kry`, `cvDiurnal_kry_bp`, `cvKrylovDemo_ls`,
  `cvKrylovDemo_prec`, `cvParticle_dns`, `cvPendulum_dns`,
  `cvRoberts_dns`, `cvRoberts_dns_constraints`, `cvRoberts_dns_negsol`,
  `cvRoberts_dns_uw`, `cvRocket_dns`, and `cvVdp_auto_nls`: serial CVODE
  examples.
- `cvAdvDiff_bndL` and `cvRoberts_dnsL`: serial LAPACK examples.
- `cvRoberts_block_klu` and `cvRoberts_klu`: serial KLU examples.
- `cvRoberts_sps`: serial SuperLU-MT example.
- `cvAdvDiff_diag_p`, `cvAdvDiff_non_p`, `cvDiurnal_kry_bbd_p`, and
  `cvDiurnal_kry_p`: MPI parallel examples.
- `cvDiurnal_kry_mpimanyvec`: MPI/MPIManyVector example.
- `cvAdvDiff_non_ph`: MPI/HYPRE example.
- `cvAdvDiff_petsc` and `cv_petsc_ex7`: PETSc examples.
- `cvAdvDiff_bnd_omp`: OpenMP example.
- `cvAdvDiff_kry_ompdev`: OpenMP device-offloading example.

Optional-backend examples are built only when the corresponding SUNDIALS
feature is enabled. `cvParticle_dns` and `cvPendulum_dns` include plotting
scripts, while `cvRoberts_dns` includes its statistics data file.
