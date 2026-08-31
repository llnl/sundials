.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added the :cmakeop:`SUNDIALS_ENABLE_DEPRECATED_WARNINGS` CMake option to allow
users to disable compiler warnings emitted by deprecated SUNDIALS API
annotations.

Added the function :c:func:`N_VSetDeviceArrayPointer` to the N_Vector API to set
the device data pointer for vector implementations with a device memory space.

Added Butcher tables for the classical RK4
(:c:enumerator:`ARKODE_KUTTA_RK4a_4_4`) and 3/8-rule
(:c:enumerator:`ARKODE_KUTTA_RK4b_4_4`) methods.

**Bug Fixes**

Fixed bug in SUNNonlinearSolver_Auto which resulted in premature switch to Newton from fixed point
due to convergence rate check occurring after only one iteration.

Fixed duplicate keys in IDA and IDAS logging output from consistent initial
condition solves and order selection diagnostics.

**Deprecation Notices**
