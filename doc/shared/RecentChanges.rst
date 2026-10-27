.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Breaking Changes**

Functions marked with ``SUNDIALS_DEPRECATED_EXPORT_MSG`` or ``[[deprecated]]``
in prior releases have been removed and along with corresponding
``_deprecated.h`` header files.

Deprecated Runge--Kutta table constants in ARKODE have been removed.

Include guards in all SUNDIALS header files have been updated to use the naming
convention ``SUNDIALS_[FILENAME]``. Any codes relying on these macros will need
to update the names accordingly.

**Major Features**

**New Features and Enhancements**

Added the ``SUNDIALS_ENABLE_CODE_COVERAGE`` CMake option and ``codecov`` target for generating code
coverage reports.

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

**Deprecation Notices**
