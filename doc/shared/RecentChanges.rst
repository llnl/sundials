.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

**Bug Fixes**

Fixed a CMake bug where the SuperLU_MT interface would not be built and
installed without setting the ``SUPERLUMT_WORKS`` option to ``TRUE``.

Fixed the embedded coefficients for the ``ARKODE_TSITOURAS_7_4_5`` Butcher
table.

**Deprecation Notices**

Several CMake options have been deprecated in favor of namespaced versions
prefixed with ``SUNDIALS_`` to avoid naming collisions in applications that
include SUNDIALS directly within their CMake builds. Additionally, a consistent
naming convention (``SUNDIALS_ENABLE``) is now used for all boolean options. The
table below lists the old CMake option names and the new replacements.

+-------------------------------------------+---------------------------------------------------------+
| Old Option                                | New Option                                              |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_MPI``                            | :cmakeop:`SUNDIALS_ENABLE_MPI`                          |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_OPENMP``                         | :cmakeop:`SUNDIALS_ENABLE_OPENMP`                       |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_OPENMP_DEVICE``                  | :cmakeop:`SUNDIALS_ENABLE_OPENMP_DEVICE`                |
+-------------------------------------------+---------------------------------------------------------+
| ``OPENMP_DEVICE_WORKS``                   | :cmakeop:`SUNDIALS_ENABLE_OPENMP_DEVICE_CHECKS`         |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_PTHREAD``                        | :cmakeop:`SUNDIALS_ENABLE_PTHREAD`                      |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_CUDA``                           | :cmakeop:`SUNDIALS_ENABLE_CUDA`                         |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_HIP``                            | :cmakeop:`SUNDIALS_ENABLE_HIP`                          |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_SYCL``                           | :cmakeop:`SUNDIALS_ENABLE_SYCL`                         |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_LAPACK``                         | :cmakeop:`SUNDIALS_ENABLE_LAPACK`                       |
+-------------------------------------------+---------------------------------------------------------+
| ``LAPACK_WORKS``                          | :cmakeop:`SUNDIALS_ENABLE_LAPACK_CHECKS`                |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_GINKGO``                         | :cmakeop:`SUNDIALS_ENABLE_GINKGO`                       |
+-------------------------------------------+---------------------------------------------------------+
| ``GINKGO_WORKS``                          | :cmakeop:`SUNDIALS_ENABLE_GINKGO_CHECKS`                |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_MAGMA``                          | :cmakeop:`SUNDIALS_ENABLE_MAGMA`                        |
+-------------------------------------------+---------------------------------------------------------+
| ``MAGMA_WORKS``                           | :cmakeop:`SUNDIALS_ENABLE_MAGMA_CHECKS`                 |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_SUPERLUDIST``                    | :cmakeop:`SUNDIALS_ENABLE_SUPERLUDIST`                  |
+-------------------------------------------+---------------------------------------------------------+
| ``SUPERLUDIST_WORKS``                     | :cmakeop:`SUNDIALS_ENABLE_SUPERLUDIST_CHECKS`           |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_SUPERLUMT``                      | :cmakeop:`SUNDIALS_ENABLE_SUPERLUMT`                    |
+-------------------------------------------+---------------------------------------------------------+
| ``SUPERLUMT_WORKS``                       | :cmakeop:`SUNDIALS_ENABLE_SUPERLUMT_CHECKS`             |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_KLU``                            | :cmakeop:`SUNDIALS_ENABLE_KLU`                          |
+-------------------------------------------+---------------------------------------------------------+
| ``KLU_WORKS``                             | :cmakeop:`SUNDIALS_ENABLE_KLU_CHECKS`                   |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_HYPRE``                          | :cmakeop:`SUNDIALS_ENABLE_HYPRE`                        |
+-------------------------------------------+---------------------------------------------------------+
| ``HYPRE_WORKS``                           | :cmakeop:`SUNDIALS_ENABLE_HYPRE_CHECKS`                 |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_PETSC``                          | :cmakeop:`SUNDIALS_ENABLE_PETSC`                        |
+-------------------------------------------+---------------------------------------------------------+
| ``PETSC_WORKS``                           | :cmakeop:`SUNDIALS_ENABLE_PETSC_CHECKS`                 |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_TRILINOS``                       | :cmakeop:`SUNDIALS_ENABLE_TRILINOS`                     |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_RAJA``                           | :cmakeop:`SUNDIALS_ENABLE_RAJA`                         |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_XBRAID``                         | :cmakeop:`SUNDIALS_ENABLE_XBRAID`                       |
+-------------------------------------------+---------------------------------------------------------+
| ``XBRAID_WORKS``                          | :cmakeop:`SUNDIALS_ENABLE_XBRAID_CHECKS`                |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_ONEMKL``                         | :cmakeop:`SUNDIALS_ENABLE_ONEMKL`                       |
+-------------------------------------------+---------------------------------------------------------+
| ``ONEMKL_WORKS``                          | :cmakeop:`SUNDIALS_ENABLE_ONEMKL_CHECKS`                |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_CALIPER``                        | :cmakeop:`SUNDIALS_ENABLE_CALIPER`                      |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_ADIAK``                          | :cmakeop:`SUNDIALS_ENABLE_ADIAK`                        |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_KOKKOS``                         | :cmakeop:`SUNDIALS_ENABLE_KOKKOS`                       |
+-------------------------------------------+---------------------------------------------------------+
| ``KOKKOS_WORKS``                          | :cmakeop:`SUNDIALS_ENABLE_KOKKOS_CHECKS`                |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_KOKKOS_KERNELS``                 | :cmakeop:`SUNDIALS_ENABLE_KOKKOS_KERNELS`               |
+-------------------------------------------+---------------------------------------------------------+
| ``KOKKOS_KERNELS_WORKS``                  | :cmakeop:`SUNDIALS_ENABLE_KOKKOS_KERNELS_CHECKS`        |
+-------------------------------------------+---------------------------------------------------------+
| ``BUILD_BENCHMARKS``                      | :cmakeop:`SUNDIALS_ENABLE_BENCHMARKS`                   |
+-------------------------------------------+---------------------------------------------------------+
| ``BENCHMARKS_INSTALL_PATH``               | :cmakeop:`SUNDIALS_BENCHMARKS_INSTALL_PATH`             |
+-------------------------------------------+---------------------------------------------------------+
| ``SUNDIALS_BENCHMARK_OUTPUT_DIR``         | :cmakeop:`SUNDIALS_BENCHMARKS_OUTPUT_DIR`               |
+-------------------------------------------+---------------------------------------------------------+
| ``SUNDIALS_BENCHMARK_CALIPER_OUTPUT_DIR`` | :cmakeop:`SUNDIALS_BENCHMARKS_CALIPER_OUTPUT_DIR`       |
+-------------------------------------------+---------------------------------------------------------+
| ``SUNDIALS_BENCHMARK_NUM_CPUS``           | :cmakeop:`SUNDIALS_BENCHMARKS_NUM_CPUS`                 |
+-------------------------------------------+---------------------------------------------------------+
| ``SUNDIALS_BENCHMARK_NUM_GPUS``           | :cmakeop:`SUNDIALS_BENCHMARKS_NUM_GPUS`                 |
+-------------------------------------------+---------------------------------------------------------+
