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

Added :c:func:`SUNAdjointCheckpointScheme_Enable` to temporarily disable and
re-enable checkpointing while preserving the configured checkpointing strategy.

Changed the :c:func:`SUNAdjointStepper_ReInit` signature to remove the forward
initial time and state arguments and add the final forward step index. This is
a breaking API change; callers must now reinitialize and evolve the forward
stepper before calling :c:func:`SUNAdjointStepper_ReInit`.

**Bug Fixes**

**Deprecation Notices**
