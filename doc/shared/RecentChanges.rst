.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added CUDA support to the sundials4py Python interface. CUDA N_Vectors can be
used with CUDA device arrays from CuPy, PyTorch and JAX when built with
:cmakeop:`SUNDIALS_ENABLE_CUDA` set to ``ON``.

sundials4py now has a ``N_VGetNumpyArray`` function which should be used instead
of ``N_VGetArrayPointer``.

Added the :cmakeop:`SUNDIALS_ENABLE_DEPRECATED_WARNINGS` CMake option to allow
users to disable compiler warnings emitted by deprecated SUNDIALS API
annotations.

Added the function :c:func:`N_VSetDeviceArrayPointer` to the N_Vector API to set
the device data pointer for vector implementations with a device memory space.

The SUNLinearSolver interface to SuperLU_DIST now supports single precision.

Added Butcher tables for the classical RK4
(:c:enumerator:`ARKODE_KUTTA_RK4a_4_4`) and 3/8-rule
(:c:enumerator:`ARKODE_KUTTA_RK4b_4_4`) methods.

The KLU SUNLinearSolver is now available in sundials4py.

**Bug Fixes**

Fixed a bug where the factor provided by ``ARKodeSetEpsLin`` was scaled by 0.1.
To restore the original behavior, call ``ARKodeSetEpsLin`` with an argument of
0.005.

Fixed bug in SUNNonlinearSolver_Auto which resulted in premature switch to Newton from fixed point
due to convergence rate check occurring after only one iteration.

Fixed duplicate keys in IDA and IDAS logging output from consistent initial
condition solves and order selection diagnostics.

**Deprecation Notices**
