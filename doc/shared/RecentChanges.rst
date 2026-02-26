.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

ARKODE now allows users to supply functions that will be called before each internal
time step, after each successful time step, after each failed time step, before
right-hand side routines are called on an updated state, and/or once each internal
stage is computed (:c:func:`ARKodeSetPreprocessStepFn`,
:c:func:`ARKodeSetPostprocessStepFn`, :c:func:`ARKodeSetPostprocessStepFailFn`,
:c:func:`ARKodeSetPreRHSProcessFn`, and :c:func:`ARKodeSetPostprocessStageFn`).
These are considered **advanced** functions, as they should treat the state vector as
read-only, otherwise all theoretical guarantees of solution accuracy and stability
will be lost.

The default numbers of stages for the SSP Runge--Kutta methods :c:enumerator:`ARKODE_LSRK_SSP_S_2`
and :c:enumerator:`ARKODE_LSRK_SSP_S_3` in LSRKStep were changed from 10 and 9, respectively, to
their minimum allowable values of 2 and 4.  Users may revert to the previous values by calling
:c:func:`LSRKStepSetNumSSPStages`.

**Bug Fixes**

Fixed a CMake bug where the SuperLU_MT interface would not be built and
installed without setting the ``SUPERLUMT_WORKS`` option to ``TRUE``.

Fixed a bug in logging output from ARKODE, where for some time stepping modules, the
the current "time" output in the logger was incorrect.

**Deprecation Notices**
