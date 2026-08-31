.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added the :cmakeop:`SUNDIALS_ENABLE_DEPRECATED_WARNINGS` CMake option to allow
users to disable compiler warnings emitted by deprecated SUNDIALS API
annotations.

Added the function :c:func:`N_VSetDeviceArrayPointer` to the N_Vector API to set
the device data pointer for vector implementations with a device memory space.

Added support for ARKODE, CVODE, and CVODES to compute internal
difference-quotient Jacobian approximations with user-provided CSC sparse
matrix patterns.

**Bug Fixes**

Fixed bug in SUNNonlinearSolver_Auto which resulted in premature switch to Newton from fixed point
due to convergence rate check occurring after only one iteration.

**Deprecation Notices**
