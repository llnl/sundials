.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added the function :c:func:`SUNLogger_SetQueueAndFlushMsgFns` to allow for
user-defined functions to queue and flush log messages.

Updated ``examples/cvode/petsc/cv_petsc_ex7.c`` to support PETSc 3.25.0.

**Bug Fixes**

Fixed memory leaks in CVODES, IDAS, and KINSOL in the unlikely event of a failed
``malloc``.

Fixed minor bug in reporting the maximum number of stages in
:c:func:`ARKodeGetStageIndex` when running SSP methods in LSRKStep.

Removed duplicate logging output that would cause the Python logging tools to fail
with a repeated key error.

Fixed empty ``elseif()`` cases in the CMake files for the Fortran interfaces to
the ManyVector and MPIPlusX vectors which could results in a missing include
path when compiling if an MPI compiler wrapper is not found.

**Deprecation Notices**
