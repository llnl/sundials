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
close. Previously, IDA(S) would return ``IDA_ILL_INPUT``.

Fixed a bug in ARKODE, CVODE(S), and IDA(S) where the linear solver counters
were not reinitialized until the next call to advance the system. As such,
non-zero linear solver statistics could be returned if retrieving or printing
linear solver counters between the initialization and the next call to advance
the system.

The interface to Ginkgo batched linear solvers has been updated to fix build
errors when using 64-bit index types. Note, only the batched dense matrix in
Ginkgo is currently compatible with 64-bit indexing (as of Ginkgo 1.10).

The SPRKStep module now accounts for zero coefficients in the SPRK tables,
eliminating extraneous function evaluations.

A bug preventing a user supplied :c:func:`SUNStepper_ResetCheckpointIndex`
function from being called was fixed.

The Kokkos N_Vector now properly handles unmanaged views. Previously, if a
Kokkos ``N_Vector`` was created from an unmanaged view, the view would become a
managed view and the data would be freed unexpectedly.

**Deprecation Notices**
