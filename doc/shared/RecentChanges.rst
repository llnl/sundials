.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Updated the MRIStep module in ARKODE to support *extended
Super Time Stepping* (ExtSTS) methods for multi-physics simulations involving
parabolic components, through the :c:func:`MRIStepExtSTSCreate` function.
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

**Bug Fixes**

Fixed bug in SUNNonlinearSolver_Auto which resulted in premature switch to Newton from fixed point
due to convergence rate check occurring after only one iteration.

**Deprecation Notices**
