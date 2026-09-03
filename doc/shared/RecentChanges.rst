.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**New Features and Enhancements**

Updated the MRIStep module in ARKODE to support *extended Super Time Stepping*
(ExtSTS) methods for multi-physics simulations involving parabolic components,
through the :c:func:`MRIStepCreateExtSTS` function.  See the documentation
section on :ref:`the mathematical methods in MRIStep
<ARKODE.Mathematics.MRIStep.ExtSTS>` for more details on the structure of ExtSTS
methods, and the documentation section on :ref:`a skeleton of usage for MRIStep
<ARKODE.Usage.MRIStep.Skeleton-ExtSTS>` for details on its usage.

Added Butcher tables for the classical RK4
(:c:enumerator:`ARKODE_KUTTA_RK4a_4_4`) and 3/8-rule
(:c:enumerator:`ARKODE_KUTTA_RK4b_4_4`) methods.

Added CUDA support to the sundials4py Python interface. CUDA N_Vectors can be
used with CUDA device arrays from CuPy, PyTorch and JAX when built with
:cmakeop:`SUNDIALS_ENABLE_CUDA` set to ``ON``.

The KLU SUNLinearSolver is now available in sundials4py.

sundials4py now has a ``N_VGetNumpyArray`` function which should be used instead
of ``N_VGetArrayPointer``.

The SUNLinearSolver interface to SuperLU_DIST now supports single precision.

Added the function :c:func:`N_VSetDeviceArrayPointer` to the N_Vector API to set
the device data pointer for vector implementations with a device memory space.

Added the :cmakeop:`SUNDIALS_ENABLE_DEPRECATED_WARNINGS` CMake option to allow
users to disable compiler warnings emitted by deprecated SUNDIALS API
annotations.

**Bug Fixes**

Fixed a bug where the factor provided by ``ARKodeSetEpsLin`` was scaled by 0.1.
To restore the original behavior, call ``ARKodeSetEpsLin`` with an argument of
0.005.

Fixed bug in SUNNonlinearSolver_Auto which resulted in premature switch to
Newton from fixed point due to convergence rate check occurring after only one
iteration.

Fixed duplicate keys in IDA and IDAS logging output from consistent initial
condition solves and order selection diagnostics.

**Deprecation Notices**

Renamed the ``ARKODE_ARK2_ERK_3_1_2``, ``ARKODE_ARK2_DIRK_3_1_2``,
``ARKODE_ASCHER_ERK_3_1_2``, and ``ARKODE_ASCHER_SDIRK_3_1_2`` Butcher tables to
:c:enumerator:`ARKODE_GKC21_ERK_3_1_2`,
:c:enumerator:`ARKODE_GKC21_ESDIRK_3_1_2`,
:c:enumerator:`ARKODE_ARS222_ERK_3_1_2`, and
:c:enumerator:`ARKODE_ARS222_ESDIRK_3_1_2`, respectively, to reflect the
original inventors of each method.  The previous table names are deprecated and
will be removed in a future release.

Renamed the ``ARKODE_IMEX_MRI_GARK_ARK2`` and
``ARKODE_IMEX_MRI_GARK_ASCHER_ARK2`` multirate coupling tables to
``ARKODE_IMEX_MRI_GARK_GKC21`` and ``ARKODE_IMEX_MRI_GARK_ARS222`` to reflect
the original inventors of the base Runge--Kutta tables on which these are based.
The previous table names are deprecated and will be removed in a future release.
