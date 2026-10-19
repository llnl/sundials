.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Breaking Changes**

Functions marked with ``SUNDIALS_DEPRECATED_EXPORT_MSG`` or ``[[deprecated]]``
in prior releases have been removed and along with corresponding
``_deprecated.h`` header files.

Deprecated Runge--Kutta table constants in ARKODE have been removed.

Include guards in all SUNDIALS header files have been updated to use the naming
convention ``SUNDIALS_[FILENAME]``. Any codes relying on these macros will need
to update the names accordingly.

**Major Features**

**New Features and Enhancements**

Reorganized the ARKODE, CVODE(S), IDA(S), and KINSOL examples by problem with
subdirectories for each problem variant. The subdirectories contain the example
source, input files, scripts, and reference output. The old example Makefiles
are no longer installed; use the installed CMakeLists.txt files instead.

The header-only Kokkos and Ginkgo N_Vector, SUNMatrix, and SUNLinearSolver
targets are now exported by installed SUNDIALS CMake packages. The installed
Ginkgo targets also propagate their Ginkgo dependency to consumers.

Added the utility function, :c:func:`SUNFileFlush` for flushing file
pointers. This is useful when using the Fortran 2003 interfaces.

Updated :c:func:`ARKodeGetReturnFlagName`, :c:func:`CVodeGetReturnFlagName`,
:c:func:`CVDiagGetReturnFlagName`, :c:func:`IDAGetReturnFlagName`, and
:c:func:`KINGetReturnFlagName` to return ``const char*`` rather than ``char*``.
The returned string should no longer be freed by the user.

Updated the shared embedding coefficients for the
``ARKODE_SSP_LSPUM_ERK_3_1_2`` and ``ARKODE_SSP_LSPUM_SDIRK_3_1_2`` methods to
improve their temporal accuracy estimates.  Users who wish to revert to the
previous coefficients can load the current table(s) using
:c:func:`ARKodeButcherTable_LoadERK` or :c:func:`ARKodeButcherTable_LoadDIRK`
and replace the internal ``d`` coefficients in the returned table with the
previous values, [ 4247.0/9840.0, 21343.0/108240.0, 49.0/132.0 ].

**Bug Fixes**

Fixed a segfault that could occur when handling errors without a ``SUNContext``
provided.

Fixed a bug in ``FindMAGMA.cmake`` which didn't allow use of MAGMA versions with
multiple digits in an identifier.

Fixed a bug in the LSRKStep module where the step post-processing function for
the RKL step was called with an incorrect time argument.

Fixed a bug in :c:func:`LSRKStepSetDomEigFrequency` where a negative ``nsteps``
value was stored as the update frequency instead of restoring the default value.

**Deprecation Notices**

Public declarations of ``SUNAdaptController``, ``N_Vector``,
``SUNLinearSolver``, and ``SUNNonlinearSolver`` operations specific to a
particular implementation are deprecated and have been moved into new headers
with ``*_deprecated.h`` suffixes. Corresponding generic operations should be
used instead and will be required in version 8.1.0 when the deprecated headers
are removed. For example, ``N_VScale_Serial`` is now deprecated and should be
replaced with ``N_VScale``. To continue using ``N_VScale_Serial``, include
``include/nvector/nvector_serial_deprecated.h``. A full table of the deprecated
operations and their new headers is provided below.

+------------------------+----------------------------------------------+-----------------------------------------------------------------------------------------------------+
| Component              | Deprecated function families                 | New headers to restore deprecated operations                                                        |
+========================+==============================================+=====================================================================================================+
| ``SUNAdaptController`` | 21 ``SUNAdaptController_*_<implementation>`` | ``sunadaptcontroller_imexgus_deprecated.h``, ``sunadaptcontroller_mrihtol_deprecated.h``,           |
|                        | operations                                   | ``sunadaptcontroller_soderlind_deprecated.h``                                                       |
+------------------------+----------------------------------------------+-----------------------------------------------------------------------------------------------------+
| ``N_Vector``           | 628 ``N_V*_<implementation>`` operations     | ``nvector_cuda_deprecated.h``, ``nvector_hip_deprecated.h``, ``nvector_manyvector_deprecated.h``,   |
|                        |                                              | ``nvector_mpimanyvector_deprecated.h``, ``nvector_mpiplusx_deprecated.h``,                          |
|                        |                                              | ``nvector_openmp_deprecated.h``, ``nvector_openmpdev_deprecated.h``,                                |
|                        |                                              | ``nvector_parallel_deprecated.h``, ``nvector_parhyp_deprecated.h``, ``nvector_petsc_deprecated.h``, |
|                        |                                              | ``nvector_pthreads_deprecated.h``, ``nvector_raja_deprecated.h``, ``nvector_serial_deprecated.h``,  |
|                        |                                              | ``nvector_sycl_deprecated.h``, ``nvector_trilinos_deprecated.h``                                    |
+------------------------+----------------------------------------------+-----------------------------------------------------------------------------------------------------+
| ``SUNLinearSolver``    | 140 ``SUNLinSol*_<implementation>``          | ``sunlinsol_band_deprecated.h``, ``sunlinsol_cusolversp_batchqr_deprecated.h``,                     |
|                        | operations                                   | ``sunlinsol_dense_deprecated.h``, ``sunlinsol_klu_deprecated.h``,                                   |
|                        |                                              | ``sunlinsol_lapackband_deprecated.h``, ``sunlinsol_lapackdense_deprecated.h``,                      |
|                        |                                              | ``sunlinsol_magmadense_deprecated.h``, ``sunlinsol_onemkldense_deprecated.h``,                      |
|                        |                                              | ``sunlinsol_pcg_deprecated.h``, ``sunlinsol_spbcgs_deprecated.h``,                                  |
|                        |                                              | ``sunlinsol_spfgmr_deprecated.h``, ``sunlinsol_spgmr_deprecated.h``,                                |
|                        |                                              | ``sunlinsol_sptfqmr_deprecated.h``, ``sunlinsol_superludist_deprecated.h``,                         |
|                        |                                              | ``sunlinsol_superlumt_deprecated.h``                                                                |
+------------------------+----------------------------------------------+-----------------------------------------------------------------------------------------------------+
| ``SUNNonlinearSolver`` | 50 ``SUNNonlinSol*_<implementation>``        | ``sunnonlinsol_auto_deprecated.h``, ``sunnonlinsol_fixedpoint_deprecated.h``,                       |
|                        | operations                                   | ``sunnonlinsol_newton_deprecated.h``, ``sunnonlinsol_petscsnes_deprecated.h``                       |
+------------------------+----------------------------------------------+-----------------------------------------------------------------------------------------------------+
