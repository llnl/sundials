.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**New Features and Enhancements**

Add `__float128` support with quadmath dependency and ostream integration.

:c:func:`ARKodeSetCFLFraction` now allows ``cfl_frac`` to be greater than or
equal to one.

Added an option to enable compensated summation of the time accumulator for all
of ARKODE. This was previously only an option for the SPRKStep module. The new
function to call to enable this is :c:func:`ARKodeSetUseCompensatedSums`.

**Bug Fixes**

Fixed segfaults in :c:func:`CVodeAdjInit` and :c:func:`IDAAdjInit` when called
after adjoint memory has been freed.

Fixed a CMake bug that would cause the Caliper compile test to fail at configure
time.

Fixed a bug in the CVODE/CVODES :c:func:`CVodeSetEtaFixedStepBounds` function
which disallowed setting ``eta_min_fx`` or ``eta_min_fx`` to 1.

:c:func:`SUNAdjointStepper_PrintAllStats` was reporting the wrong quantity for
the number of "recompute passes" and has been fixed.

**Deprecation Notices**

The :c:func:`SPRKStepSetUseCompensatedSums` function has been deprecated. Use
the :c:func:`ARKodeSetUseCompensatedSums` function instead.
