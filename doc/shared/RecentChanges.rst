.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added the utility function, :c:func:`SUNFileFlush` for flushing file
pointers. This is useful when using the Fortran 2003 interfaces.

**Bug Fixes**

Fixed a bug in ``FindMAGMA.cmake`` which didn't allow use of MAGMA versions with 
multiple digits in an identifier.

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
| Component              | Deprecated function families                 | New headers to restore deprecated opertions                                                         |
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
