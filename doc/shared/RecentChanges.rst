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

Added :c:func:`MRIStepSetUseExtSTSErrorEstimate` (command-line key
``arkid.use_extsts_error_estimate``) to include the embedded error estimates of
the super-time-stepping (STS) methods in the temporal error estimate of ExtSTS
methods created with :c:func:`MRIStepCreateExtSTS`, for use when the STS error
may dominate.  The option is disabled by default.

**Bug Fixes**

Fixed a segfault that could occur when handling errors without a ``SUNContext``
provided.

Fixed a bug in ``FindMAGMA.cmake`` which didn't allow use of MAGMA versions with 
multiple digits in an identifier.

Fixed a bug in the LSRKStep module where the step post-processing function for
the RKL step was called with an incorrect time argument.

Fixed a bug in :c:func:`LSRKStepSetDomEigFrequency` where a negative ``nsteps``
value was stored as the update frequency instead of restoring the default value.

Fixed the ``ark_adr1d_extsts`` example to select an explicit MRI method by
default when run with ``--no-reaction``, and corrected its input check for
``--build_mri_table``.

**Deprecation Notices**
