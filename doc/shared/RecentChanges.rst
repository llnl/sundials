.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

Replaced the ARKODE ``MRIStepInnerStepper`` interface with :c:type:`SUNStepper`.
:c:func:`MRIStepCreate` now accepts a :c:type:`SUNStepper`, which can created by
:c:func:`ARKodeCreateSUNStepper` or by manually wrapping an IVP solving
procedure. When porting a custom ``MRIStepInnerStepper`` to a
:c:type:`SUNStepper`, most operations have the same signature, but forcing is
handled slightly differently. A :c:type:`SUNStepper` is responsible for storing
forcing data set via :c:func:`SUNStepper_SetForcing` and applying it in
:c:func:`SUNStepper_Evolve`, possibly with the help of the new
:c:func:`SUNStepper_AddForcing` function. See
``examples/arkode/CXX_serial/ark_advection_diffusion_reaction.cpp`` for an
example. :c:func:`SUNStepper_Evolve` and :c:type:`SUNStepperEvolveFn` now return
an integer status where zero is success, positive values are recoverable
failures, and negative values are fatal failures. Added accumulated error
get/reset and relative tolerance operations to :c:type:`SUNStepper` to support
:c:type:`SUNAdaptController_MRIHTol`.

Added discrete adjoint support for explicit MRI-GARK methods in MRIStep with the
:c:func:`MRIStepCreateAdjointStepper` function which enables sensitivity
analysis capabilities for ODEs with disparate time scales. This uses a new
:c:type:`MRIStepInnerAdjointProblem` object which augments the terminal state,
right-hand side function, and user data of the inner adjoint integrator to
include terms need to account for the forcing MRIStep applies during the forward
integration.

**New Features and Enhancements**

Added :c:func:`SUNAdjointCheckpointScheme_Enable` to temporarily disable and
re-enable checkpointing while preserving the configured checkpointing strategy.

Changed the :c:func:`SUNAdjointStepper_ReInit` signature to remove the forward
initial time and state arguments and add the final forward step index. This is
a breaking API change; callers must now reinitialize and evolve the forward
stepper before calling :c:func:`SUNAdjointStepper_ReInit`.

Changed the argument of :c:type:`SUNAdjointCheckpointSchemeDestroyFn` from
``SUNAdjointCheckpointScheme*`` to ``SUNAdjointCheckpointScheme``.
Implementations of this function only need to free content associated with the
derived c:type:`SUNAdjointCheckpointScheme`.

Added the utility function, :c:func:`SUNFileFlush` for flushing file
pointers. This is useful when using the Fortran 2003 interfaces.

Updated :c:func:`ARKodeGetReturnFlagName`, :c:func:`CVodeGetReturnFlagName`,
:c:func:`CVDiagGetReturnFlagName`, :c:func:`IDAGetReturnFlagName`, and
:c:func:`KINGetReturnFlagName` to return ``const char*`` rather than ``char*``.
The returned string should no longer be freed by the user.

**Bug Fixes**

Fixed a bug in the ERKStep and ARKStep discrete adjoints where it would take an
extra backward step when the last forward step was slightly truncated to align
with the end time.

Fixed a segfault that could occur when handling errors without a ``SUNContext``
provided.

Fixed a bug in ``FindMAGMA.cmake`` which didn't allow use of MAGMA versions with 
multiple digits in an identifier.

Fixed a bug in the LSRKStep module where the step post-processing function for
the RKL step was called with an incorrect time argument.

**Deprecation Notices**
