.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added functions to CVODE(S) and IDA(S) to set the maximum number of inequality
constraint failures in a step attempt (:c:func:`CVodeSetMaxNumConstraintFails`
and :c:func:`IDASetMaxNumConstraintFails`) and to retrieve the total number of
failed step attempts due to an inequality constraint violation
(:c:func:`CVodeGetNumConstraintFails` and
:c:func:`IDAGetNumConstraintFails`). As a result, constraint failures are no
longer included in the number of step failures due to a solver failure (i.e.,
the values returned by :c:func:`CVodeGetNumStepSolveFails` and
:c:func:`IDAGetNumStepSolveFails`).

**Bug Fixes**

Fixed a bug in the CVODE(S) inequality constraint handling where the predicted
state was used to compute the step size reduction factor which could lead to an
insufficient reduction in the step size or, when the prediction violates the
constraints, an infinitely large step size in the next step attempt.

The SPRKStep module now accounts for zero coefficients in the SPRK tables,
eliminating extraneous function evaluations.

A bug preventing a user supplied :c:func:`SUNStepper_ResetCheckpointIndex`
function from being called was fixed.

In CVODES and IDA, added missing return flag names to
:c:func:`CVodeGetReturnFlagName` and :c:func:`IDAGetReturnFlagName`,
respectively.

**Deprecation Notices**
