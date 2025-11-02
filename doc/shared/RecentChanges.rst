.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

**Bug Fixes**

On the initial time step, ARKODE and CVODE(S) will now return ``ARK_TOO_CLOSE``
or ``CV_TOO_COLSE``, respectively, when the requested output time is the same as
the initial time (or within numerical roundoff of the initial time) with a
user-supplied initial step size. Before a ``TOO_CLOSE`` error would only be
returned when internally estimating the initial step size. In IDA(S), added a
``IDA_TOO_CLOSE`` return value for when the initial and output time are too
close. Previosuly, IDA(S) would return ``IDA_ILL_INPUT``.

Fixed a bug in ARKODE, CVODE(S), and IDA(S) where the linear solver counters
were not reinitialized until the next call to advance the system. As such,
non-zero linear solver statistics could be returned if retrieving or printing
linear solver counters between the initialization and the next call to advance
the system.

The SPRKStep module now accounts for zero coefficients in the SPRK tables,
eliminating extraneous function evaluations.

A bug preventing a user supplied :c:func:`SUNStepper_ResetCheckpointIndex`
function from being called was fixed.

**Deprecation Notices**
