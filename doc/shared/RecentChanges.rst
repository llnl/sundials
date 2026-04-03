.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**New Features and Enhancements**

The default number of stages for the SSP Runge-Kutta methods
:c:enumerator:`ARKODE_LSRK_SSP_S_2` and :c:enumerator:`ARKODE_LSRK_SSP_S_3` in
LSRKStep were changed from 10 and 9, respectively, to their minimum allowable
values of 2 and 4. Users may revert to the previous values by calling
:c:func:`LSRKStepSetNumSSPStages`.

Added the optional function :c:func:`ARKodeInit` to ARKODE to enable data
allocation before the first call to :c:func:`ARKodeEvolve` (but after all other
optional input routines have been called), to support users who measure memory
usage before beginning a simulation.

Added the function :c:func:`ARKodeGetStageIndex` that returns the index of the
stage currently being processed, and the total number of stages in the method,
for users who wish to compute auxiliary quantities in their IVP right-hand side
functions during some stages and not others (e.g., in all but the first or last
stage).

Added the functions :c:func:`ARKodeGetLastTime` and :c:func:`ARKodeGetLastState`
to return the last successful time and state achieved by ARKODE, respectively.

ARKODE now allows users to supply functions that will be called before each
internal time step attempt (:c:func:`ARKodeSetPreStepFn`), after each successful
time step (:c:func:`ARKodeSetPostStepFn`), before right-hand side routines are
called on an updated state (:c:func:`ARKodeSetPreRhsFn`), and/or once each
internal step/stage is computed (:c:func:`ARKodeSetPostprocessStepFn`/
:c:func:`ARKodeSetPostprocessStageFn`). These are considered **advanced**
functions, as they should treat the state vector as read-only, otherwise all
theoretical guarantees of solution accuracy and stability will be lost.  As a
result of these new functions, the values of multiple ARKODE return codes (e.g.,
``ARK_INTERP_FAIL``) have been updated; users who key off of the named constants
will not be affected, but users who rely on the values themselves should update
their codes accordingly.

Note to users utilizing the previously undocumented
:c:func:`ARKodeSetPostprocessStepFn` function, the supplied function is now
called on the newly computed state vector for all step attempts not just
successful steps. To obtain the previous behavior of only calling a function on
successful steps, switch to using :c:func:`ARKodeSetPostStepFn`.

Added ``SUNLogger_Set{Error,Warning,Info,Debug}File`` functions to allow setting
logger output streams with a ``FILE*``.

Updated the Kokkos N_Vector to support Kokkos 5.x versions.

**Bug Fixes**

Fixed a CMake bug where the SuperLU_MT interface would not be built and
installed without setting the ``SUPERLUMT_WORKS`` option to ``TRUE``.

Fixed the embedded coefficients for the ``ARKODE_TSITOURAS_7_4_5`` Butcher
table.

Fixed a bug in LSRKStep where an incorrect state vector could be passed to a
user-supplied dominant eigenvalue function on the first step unless the output
vector passed to :c:func:`ARKodeEvolve` contained the initial condition and when
an eigenvalue estimate is requested on the first step in a subsequent call to
:c:func:`ARKodeEvolve` unless the output vector passed contained the most
recently returned solution.

Fixed a potential bug in LSRKStep's :c:enumerator:`ARKODE_LSRK_SSP_S_3` method,
where a real number was used instead of an integer, potentially resulting in a
rounding error.

Fixed a bug in MRIStep for estimating the first "slow" time step in an adaptive
multirate calculation.

Fixed a bug in MRIStep when using a custom inner integrator that relies on the
input state being the initial condition for the fast integration rather than
retaining the result from the last inner integration or most recent reset call
and the output vector passed to :c:func:`ARKodeEvolve` does not contain the
initial condition on the first call or the last returned solution on subsequent
calls.

Added a missing call to :c:func:`SUNNonlinSolSetup` in MRIStep when using an
IMEX-MRI-SR method.

Fixed a bug in the ARKODE discrete adjoint checkpointing where an incorrect
state would be stored on the first step if the output vector passed to
:c:func:`ARKodeEvolve` did not contain the initial condition on the first call.

Removed extraneous copy of output vector when using ARKODE in ``ARK_ONE_STEP``
mode.

Removed an extraneous copy of the output vector in each step with SplittingStep.

Fixed a bug in logging output from ARKODE, where for some time stepping modules,
the current "time" output in the logger was incorrect.

Fixed a bug where passing an empty string to
``SUNLogger_Set{Error,Warning,Info,Debug}Filename`` did not disable the
corresponding logging stream `Issue #844
<https://github.com/llnl/sundials/issues/844>`__.

**Deprecation Notices**

The ``CVodeSetMonitorFn`` and ``CVodeSetMonitorFrequency`` functions have been
deprecated and will be removed in the next major release.

Several CMake options have been deprecated in favor of namespaced versions
prefixed with ``SUNDIALS_`` to avoid naming collisions in applications that
include SUNDIALS directly within their CMake builds. Additionally, a consistent
naming convention (``SUNDIALS_ENABLE``) is now used for all boolean options. The
table below lists the old CMake option names and the new replacements.

+-------------------------------------------+---------------------------------------------------------+
| Old Option                                | New Option                                              |
+-------------------------------------------+---------------------------------------------------------+
| ``BUILD_ARKODE``                          | :cmakeop:`SUNDIALS_ENABLE_ARKODE`                       |
+-------------------------------------------+---------------------------------------------------------+
| ``BUILD_CVODE``                           | :cmakeop:`SUNDIALS_ENABLE_CVODE`                        |
+-------------------------------------------+---------------------------------------------------------+
| ``BUILD_CVODES``                          | :cmakeop:`SUNDIALS_ENABLE_CVODES`                       |
+-------------------------------------------+---------------------------------------------------------+
| ``BUILD_IDA``                             | :cmakeop:`SUNDIALS_ENABLE_IDA`                          |
+-------------------------------------------+---------------------------------------------------------+
| ``BUILD_IDAS``                            | :cmakeop:`SUNDIALS_ENABLE_IDAS`                         |
+-------------------------------------------+---------------------------------------------------------+
| ``BUILD_KINSOL``                          | :cmakeop:`SUNDIALS_ENABLE_KINSOL`                       |
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
