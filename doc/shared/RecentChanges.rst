.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

Replaced the ARKODE ``MRIStepInnerStepper`` interface with :c:type:`SUNStepper`.
:c:func:`MRIStepCreate` now accepts a :c:type:`SUNStepper`, which can created by
:c:func:`ARKodeCreateSUNStepper` or manually. :c:func:`SUNStepper_Evolve` and
:c:type:`SUNStepperEvolveFn` now return an integer status where zero is success,
positive values are recoverable failures, and negative values are fatal
failures. Added accumulated error get/reset and relative tolerance operations
to `SUNStepper` to support `SUNAdaptController_MRIHTol`.

**New Features and Enhancements**

Added the :cmakeop:`SUNDIALS_ENABLE_DEPRECATED_WARNINGS` CMake option to allow
users to disable compiler warnings emitted by deprecated SUNDIALS API
annotations.

Added the function :c:func:`N_VSetDeviceArrayPointer` to the N_Vector API to set
the device data pointer for vector implementations with a device memory space.

**Bug Fixes**

**Deprecation Notices**
