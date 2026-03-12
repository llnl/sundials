.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Updated the Kokkos N_Vector to support Kokkos 5.x versions.

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
| ``BUILD_FORTRAN_MODULE_INTERFACE``        | :cmakeop:`SUNDIALS_ENABLE_FORTRAN`                      |
+-------------------------------------------+---------------------------------------------------------+
| ``SUNDIALS_BUILD_WITH_PROFILING``         | :cmakeop:`SUNDIALS_ENABLE_PROFILING`                    |
+-------------------------------------------+---------------------------------------------------------+
| ``SUNDIALS_BUILD_WITH_MONITORING``        | :cmakeop:`SUNDIALS_ENABLE_MONITORING`                   |
+-------------------------------------------+---------------------------------------------------------+
| ``SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS``  | :cmakeop:`SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS`        |
+-------------------------------------------+---------------------------------------------------------+
| ``EXAMPLES_ENABLE_C``                     | :cmakeop:`SUNDIALS_ENABLE_C_EXAMPLES`                   |
+-------------------------------------------+---------------------------------------------------------+
| ``EXAMPLES_ENABLE_CXX``                   | :cmakeop:`SUNDIALS_ENABLE_CXX_EXAMPLES`                 |
+-------------------------------------------+---------------------------------------------------------+
| ``EXAMPLES_ENABLE_F2003``                 | :cmakeop:`SUNDIALS_ENABLE_FORTRAN_EXAMPLES`             |
+-------------------------------------------+---------------------------------------------------------+
| ``EXAMPLES_ENABLE_CUDA``                  | :cmakeop:`SUNDIALS_ENABLE_CUDA_EXAMPLES`                |
+-------------------------------------------+---------------------------------------------------------+
| ``EXAMPLES_INSTALL``                      | :cmakeop:`SUNDIALS_ENABLE_EXAMPLES_INSTALL`             |
+-------------------------------------------+---------------------------------------------------------+
| ``EXAMPLES_INSTALL_PATH``                 | :cmakeop:`SUNDIALS_EXAMPLES_INSTALL_PATH`               |
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
| ``ENABLE_ALL_WARNINGS``                   | :cmakeop:`SUNDIALS_ENABLE_ALL_WARNINGS`                 |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_WARNINGS_AS_ERRORS``             | :cmakeop:`CMAKE_COMPILE_WARNING_AS_ERROR`               |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_ADDRESS_SANITIZER``              | :cmakeop:`SUNDIALS_ENABLE_ADDRESS_SANITIZER`            |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_MEMORY_SANITIZER``               | :cmakeop:`SUNDIALS_ENABLE_MEMORY_SANITIZER`             |
+-------------------------------------------+---------------------------------------------------------+
| ``ENABLE_LEAK_SANITIZER``                 | :cmakeop:`SUNDIALS_ENABLE_LEAK_SANITIZER`               |
+-------------------------------------------+---------------------------------------------------------+

Following the updated CMake options, the macros listed below have been
deprecated and replaced with versions that align with the new CMake options.

+------------------------------------------+-------------------------------------------+
| Old Macro                                | New Macro                                 |
+------------------------------------------+-------------------------------------------+
| ``SUNDIALS_BUILD_WITH_PROFILING``        | ``SUNDIALS_ENABLE_PROFILING``             |
+------------------------------------------+-------------------------------------------+
| ``SUNDIALS_BUILD_WITH_MONITORING``       | ``SUNDIALS_ENABLE_MONITORING``            |
+------------------------------------------+-------------------------------------------+
| ``SUNDIALS_BUILD_PACKAGE_FUSED_KERNELS`` | ``SUNDIALS_ENABLE_PACKAGE_FUSED_KERNELS`` |
+------------------------------------------+-------------------------------------------+
