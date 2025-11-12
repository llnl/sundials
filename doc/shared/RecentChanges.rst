.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

**Bug Fixes**

The interface to Ginkgo batched linear solvers has been updated to fix build
errors when using 64-bit index types. Note, only the batched dense matrix in
Ginkgo is currently compatible with 64-bit indexing (as of Ginkgo 1.10).

The SPRKStep module now accounts for zero coefficients in the SPRK tables,
eliminating extraneous function evaluations.

A bug preventing a user supplied :c:func:`SUNStepper_ResetCheckpointIndex`
function from being called was fixed.

**Deprecation Notices**
