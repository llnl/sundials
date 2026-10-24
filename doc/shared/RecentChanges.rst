.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added the utility function, :c:func:`SUNFileFlush` for flushing file
pointers. This is useful when using the Fortran 2003 interfaces.

Added the optional function :c:func:`SUNDomEigEstimator_SetPreRhsFn` to the
:c:type:`SUNDomEigEstimator` API to supply a right-hand side preprocessing
function that is called before each right-hand side evaluation inside the
internal difference quotient Jacobian-vector product. The function
:c:func:`SUNDomEigEstimator_SetRhsLinearizationPoint` now takes an additional
(optional) argument for the value of the right-hand side function at the
linearization point; when supplied, the estimator reuses this vector instead of
re-evaluating the right-hand side. LSRKStep now attaches its right-hand side
function to the :c:type:`SUNDomEigEstimator` through
:c:func:`SUNDomEigEstimator_SetRhs` and supplies the current right-hand side
vector, so the dominant eigenvalue estimation no longer re-evaluates the
right-hand side at the start of a step. As a result, the number of right-hand
side evaluations reported by :c:func:`LSRKStepGetNumDomEigEstRhsEvals` is now
the cumulative count returned by :c:func:`SUNDomEigEstimator_GetNumRhsEvals`.
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

Fixed a bug in the Power and Arnoldi :c:type:`SUNDomEigEstimator` implementations
where :c:func:`SUNDomEigEstimator_SetRhs` unconditionally replaced the
estimator's ``ATimes`` routine with the internal difference quotient
approximation, so a user-supplied ``ATimes`` function (e.g., one attached after
:c:func:`LSRKStepSetDomEigEstimator`) was silently discarded. The internal
approximation is now installed only when no user-supplied ``ATimes`` function is
present.
Passing ``NULL`` for the ``ATimes`` argument of
:c:func:`SUNDomEigEstimator_SetATimes` now reverts these estimators to the
internal approximation.

**Deprecation Notices**
