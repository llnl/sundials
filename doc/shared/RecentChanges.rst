.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added CUDA support to the sundials4py Python interface. CUDA N_Vectors can be
used with CUDA device arrays from CuPy, PyTorch and JAX when built with
:cmakeop:`SUNDIALS_ENABLE_CUDA` set to ``ON``.

sundials4py now has a ``N_VGetNumpyArray`` function which should be used instead
of ``N_VGetArrayPointer``.

Updated the MRIStep module in ARKODE to support *extended
Super Time Stepping* (ExtSTS) methods for multi-physics simulations involving
parabolic components, through the :c:func:`MRIStepCreateExtSTS` function.
See the documentation section on
:ref:`the mathematical methods in MRIStep <ARKODE.Mathematics.MRIStep.ExtSTS>`
for more details on the structure of ExtSTS methods, and the documentation
section on :ref:`a skeleton of usage for MRIStep <ARKODE.Usage.MRIStep.Skeleton-ExtSTS>`
for details on its usage.

Added the :cmakeop:`SUNDIALS_ENABLE_DEPRECATED_WARNINGS` CMake option to allow
users to disable compiler warnings emitted by deprecated SUNDIALS API
annotations.

Added the function :c:func:`N_VSetDeviceArrayPointer` to the N_Vector API to set
the device data pointer for vector implementations with a device memory space.

The KLU SUNLinearSolver is now available in sundials4py.

**Bug Fixes**

Fixed bug in SUNNonlinearSolver_Auto which resulted in premature switch to Newton from fixed point
due to convergence rate check occurring after only one iteration.

Fixed duplicate keys in IDA and IDAS logging output from consistent initial
condition solves and order selection diagnostics.

**Deprecation Notices**

Renamed the ``ARKODE_ARK2_ERK_3_1_2`` and ``ARKODE_ARK2_DIRK_3_1_2`` Butcher tables to
``ARKODE_GKC21_ERK_3_1_2`` and ``ARKODE_GKC21_DIRK_3_1_2`` to reflect the method's origin.
The previous table names are deprecated and will be removed in a future release.

Renamed the ``ARKODE_IMEX_MRI_GARK_ARK2`` multirate coupling table to
``ARKODE_IMEX_MRI_GARK_GKC21`` to reflect the method's origin. The previous table name is
deprecated and will be removed in a future release.
