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

Updated the shared embedding coefficients for the
``ARKODE_SSP_LSPUM_ERK_3_1_2`` and ``ARKODE_SSP_LSPUM_SDIRK_3_1_2`` methods to
improve their temporal accuracy estimates.  Users who wish to revert to the
previous coefficients can load the current table(s) using
:c:func:`ARKodeButcherTable_LoadERK` or :c:func:`ARKodeButcherTable_LoadDIRK`
and replace the internal ``d`` coefficients in the returned table with the
previous values, [ 4247.0/9840.0, 21343.0/108240.0, 49.0/132.0 ].

**Bug Fixes**

Fixed a segfault that could occur when handling errors without a ``SUNContext``
provided.

Fixed a bug in ``FindMAGMA.cmake`` which didn't allow use of MAGMA versions with
multiple digits in an identifier.

Fixed a bug in the LSRKStep module where the step post-processing function for
the RKL step was called with an incorrect time argument.

**Deprecation Notices**
