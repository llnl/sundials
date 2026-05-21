.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added the function :c:func:`SUNLogger_SetQueueAndFlushMsgFns` to allow for
user-defined functions to queue and flush log messages.

Updated ``examples/cvode/petsc/cv_petsc_ex7.c`` to support PETSc 3.25.0.

**Bug Fixes**

Fixed a minor bug where the number of required stages for STS methods 
in the LSRKStep module was incorrectly computed using the spectral 
radius instead of the real part of the Jacobian eigenvalues.

Fixed a minor bug where negative real extend of the stability region 
for the RKC method was not being properly computed, which could result 
in an underestimation of the number of stages required for stability.

Fixed a minor bug where STS methods were limited to one fewer than
the maximum allowed number of stages. STS can now use the full maximum
number of stages.

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
