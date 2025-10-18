.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added the functions :c:func:`CVodeSetMaxNumConstraintFails` and
:c:func:`CVodeGetNumConstraintFails` to set maximum number of inequality
constraint failures in a step attempt and the total number of failed step
attempts due to an inequality constraint violation, respectively.

**Bug Fixes**

Fixed a bug in the CVODE(S) inequality constraint handling where the predicted
state was used to compute the step size reduction factor which could lead to an
insufficient reduction in the step size or, when the prediction violates the
constraints, an infinitely large step size in the next step attempt.

The SPRKStep module now accounts for zero coefficients in the SPRK tables,
eliminating extraneous function evaluations.

A bug preventing a user supplied :c:func:`SUNStepper_ResetCheckpointIndex`
function from being called was fixed.

**Deprecation Notices**
