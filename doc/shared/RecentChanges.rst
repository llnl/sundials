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

sundials4py now supports implementing SUNDIALS classes in Python. The new
``CustomNVector``, ``CustomSUNMatrix``, ``CustomSUNLinearSolver``,
``CustomSUNNonlinearSolver``, ``CustomSUNDomEigEstimator``,
``CustomSUNHController``, and ``CustomSUNMRIController`` base classes may be
subclassed to provide an :c:type:`N_Vector`, :c:type:`SUNMatrix`,
:c:type:`SUNLinearSolver`, :c:type:`SUNNonlinearSolver`,
:c:type:`SUNDomEigEstimator`, or :c:type:`SUNAdaptController` implementation
written in Python, and instances of such a subclass may be passed to any
sundials4py function that expects the corresponding SUNDIALS object. See
:ref:`Python.Usage.CustomObjects` for details, and the
``cvs_custom_nonlinsol.py``, ``kin_custom_linsol.py``, and
``ark_custom_adaptcontroller.py`` examples for annotated templates.

**Bug Fixes**

Fixed IDAS staggered forward sensitivity analysis retrying a failed
converged-state residual evaluation without restoring the attempted step.
Recoverable failures now reduce the step size and count toward the existing
convergence-failure limit.

Fixed a bug in the band SUNMatrix ``SUNMatScaleAdd`` operation where entries
outside the band of the second matrix were not scaled.

Fixed :c:func:`SUNMatCopyOps` so it also copies the optional
Hermitian-transpose matrix-vector operation.

Fixed :c:func:`SUNMatClone` to return safely when an implementation's clone
operation fails and returns ``NULL``.

Fixed the sundials4py :c:func:`SUNLinSolSetPreconditioner` binding so
setup-only preconditioners are retained, the solve callback has the correct
type, and either callback may be ``None``.

Fixed sundials4py generation so
:c:func:`SUNDomEigEstimator_SetRhsLinearizationPoint` is exposed.

Fixed leaks of Python callback tables attached to native linear solvers,
nonlinear solvers, and dominant eigenvalue estimators.

Fixed the KINSOL dense difference-quotient Jacobian leaving a temporary
perturbation in the solution vector when the system function fails. The
linear-solver setup failure behavior is unchanged.

Fixed a segfault that could occur when handling errors without a ``SUNContext``
provided.

Fixed a bug in ``FindMAGMA.cmake`` which didn't allow use of MAGMA versions with
multiple digits in an identifier.

Fixed a bug in the LSRKStep module where the step post-processing function for
the RKL step was called with an incorrect time argument.

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
