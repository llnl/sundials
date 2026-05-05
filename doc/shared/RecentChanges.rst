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

Fixed memory leaks in CVODES, IDAS, and KINSOL in the unlikely event of a failed
``malloc``.

**Deprecation Notices**
