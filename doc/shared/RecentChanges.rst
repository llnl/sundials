.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Multiple minor updates were made to the ARKODE package.  We removed an extraneous
copy of the output vector when using ARKODE in ``ARK_ONE_STEP`` mode.
We added the function :c:func:`ARKodeAllocateInternalData` to ARKODE to enable
stage-related data allocation before the first call to :c:func:`ARKodeEvolve`
(but after all other optional input routines have been called), to support
users who measure memory usage before beginning a simulation.
We added the function :c:func:`ARKodeGetStageIndex` that returns the index of the
stage currently being processed, and the total number of stages in the method, for
users who must compute auxiliary quantities in their IVP right-hand side functions
during some stages and not others (e.g., in all but the first or last stage).
The default number of stages for the SSP Runge-Kutta methods :c:enumerator:`ARKODE_LSRK_SSP_S_2`
and :c:enumerator:`ARKODE_LSRK_SSP_S_3` in LSRKStep were changed from 10 and 9, respectively, to
their minimum allowable values of 2 and 4. Users may revert to the previous values by calling
:c:func:`LSRKStepSetNumSSPStages`.

ARKODE now allows users to supply functions that will be called before each internal
time step, after each successful time step, after each failed time step, before
right-hand side routines are called on an updated state, and/or once each internal
stage is computed (:c:func:`ARKodeSetPreprocessStepFn`,
:c:func:`ARKodeSetPostprocessStepFn`, :c:func:`ARKodeSetPostprocessStepFailFn`,
:c:func:`ARKodeSetPreRHSProcessFn`, and :c:func:`ARKodeSetPostprocessStageFn`).
These are considered **advanced** functions, as they should treat the state vector as
read-only, otherwise all theoretical guarantees of solution accuracy and stability
will be lost.

**Bug Fixes**

Fixed a CMake bug where the SuperLU_MT interface would not be built and
installed without setting the ``SUPERLUMT_WORKS`` option to ``TRUE``.

Fixed the embedded coefficients for the ``ARKODE_TSITOURAS_7_4_5`` Butcher
table.

Fixed a bug in logging output from ARKODE, where for some time stepping modules, the
the current "time" output in the logger was incorrect.

Fixed a potential bug in LSRKStep's :c:enumerator:`ARKODE_LSRK_SSP_S_3` method, where a real
number was used instead of an integer, potentially resulting in a rounding error.

**Deprecation Notices**
