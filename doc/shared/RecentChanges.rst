.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added the function :c:func:`SUNLogger_SetQueueAndFlushMsgFns` to allow for
user-defined functions to queue and flush log messages.

Updated ``examples/cvode/petsc/cv_petsc_ex7.c`` to support PETSc 3.25.0.

Added the function :c:func:`ARKodeSkipAdaptStopTime` to specify that
stop-time-limited steps should be disregarded when selecting step sizes
for time step adaptivity.  Added the functions
:c:func:`CVodeSkipAdaptStopTime` and :c:func:`IDASkipAdaptStopTime`
to specify that stop-time-limited steps should be disregarded when
adapting the step size and method order for CVODE(S) and IDA(S), respectively.

**Bug Fixes**

Fixed memory leaks in CVODES, IDAS, and KINSOL in the unlikely event of a failed
``malloc``.

Fixed minor bug in reporting the maximum number of stages in
:c:func:`ARKodeGetStageIndex` when running SSP methods in LSRKStep.

Removed duplicate logging output that would cause the Python logging tools to
fail with a repeated key error.

Fixed a CMake issue that prevented finding third-party libraries installed in
default search locations e.g., paths included in ``CMAKE_INSTALL_PREFIX``
(`Issue #935 <https://github.com/llnl/sundials/issues/935>`__).

Fixed empty ``elseif()`` cases in the CMake files for the Fortran interfaces to
the ManyVector and MPIPlusX vectors which could results in a missing include
path when compiling if an MPI compiler wrapper is not found.

**Deprecation Notices**
