.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Reorganized the ARKODE, CVODE(S), IDA(S), and KINSOL examples by problem with
subdirectories for each problem variant. The subdirectories contain the example
source, input files, scripts, and reference output. The old example Makefiles
are no longer installed; use the installed CMakeLists.txt files instead.

The header-only Kokkos and Ginkgo N_Vector, SUNMatrix, and SUNLinearSolver
targets are now exported by installed SUNDIALS CMake packages. The installed
Ginkgo targets also propagate their Ginkgo dependency to consumers.

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

Updated the shared embedding coefficients for the
``ARKODE_SSP_LSPUM_ERK_3_1_2`` and ``ARKODE_SSP_LSPUM_SDIRK_3_1_2`` methods to
improve their temporal accuracy estimates.  Users who wish to revert to the
previous coefficients can load the current table(s) using
:c:func:`ARKodeButcherTable_LoadERK` or :c:func:`ARKodeButcherTable_LoadDIRK`
and replace the internal ``d`` coefficients in the returned table with the
previous values, [ 4247.0/9840.0, 21343.0/108240.0, 49.0/132.0 ].

**Bug Fixes**

The Newton nonlinear solver now reports missing required callbacks during
initialization in all build configurations, preventing NULL callback crashes
across SUNDIALS packages.

Fixed IDAS staggered forward sensitivity analysis retrying a failed
converged-state residual evaluation without restoring the attempted step.
Recoverable failures now reduce the step size and count toward the existing
convergence-failure limit.

Fixed a bug in the band SUNMatrix ``SUNMatScaleAdd`` operation where entries
outside the band of the second matrix were not scaled.

Fixed the KINSOL dense difference-quotient Jacobian leaving a temporary
perturbation in the solution vector when the system function fails. The
linear-solver setup failure behavior is unchanged.

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
present. Passing ``NULL`` for the ``ATimes`` argument of
:c:func:`SUNDomEigEstimator_SetATimes` now reverts these estimators to the
internal approximation.

Fixed a bug in LSRKStep where calls
:c:func:`LSRKStepSetNumDomEigEstPreprocessIters` after the first dominant
eigenvalue estimate incorrectly set the number of preprocessing iterations to
the *initial* instead of the requested value.

Fixed a bug in :c:func:`LSRKStepSetDomEigFrequency` where a negative ``nsteps``
value was stored as the update frequency instead of restoring the default value.

Fixed the Soderlind-family and ImEx Gustafsson ``SUNAdaptController`` modules so
that a zero local error estimate no longer leads to an arbitrary step size
reduction or growth on the following steps. A zero error estimate now restarts
the controller history.

ARKODE now returns ``ARK_CONTROLLER_ERR`` when a ``SUNAdaptController`` returns
``NaN`` for the new step size or otherwise fails during a step. Previously these
failures were reported as ``ARK_ERR_FAILURE`` or ``ARK_UNRECOGNIZED_ERROR``.

**Deprecation Notices**
