.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added the utility function, :c:func:`SUNFileFlush` for flushing file
pointers. This is useful when using the Fortran 2003 interfaces.

Updated :c:func:`ARKodeGetReturnFlagName`, :c:func:`CVodeGetReturnFlagName`,
:c:func:`CVDiagGetReturnFlagName`, :c:func:`IDAGetReturnFlagName`, and
:c:func:`KINGetReturnFlagName` to return ``const char*`` rather than ``char*``.
The returned string should no longer be freed by the user.

**Bug Fixes**

Fixed a segfault that could occur when handling errors without a ``SUNContext``
provided.

Fixed a bug in ``FindMAGMA.cmake`` which didn't allow use of MAGMA versions with 
multiple digits in an identifier.

Fixed a bug in the LSRKStep module where the step post-processing function for
the RKL step was called with an incorrect time argument.

Fixed a bug in CVODES where a recoverable failure of the right-hand side
function at the converged state in the staggered sensitivity corrector
(``CV_STAGGERED`` or ``CV_STAGGERED1``) was retried without restoring the
predicted solution or reducing the step size. The next attempt then started
from the advanced time, so the step covered twice the step size and the first
half was never corrected. The failure is now handled like other recoverable
nonlinear solver failures and counts toward the maximum number of convergence
failures.

**Deprecation Notices**
