.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

SUNDIALS now has official Python interfaces! With this release, we are shipping
a **beta version** of the sundials4py Python module (created with nanobind and
litgen). sundials4py provides explicit interfaces to most features of SUNDIALS.
See the :ref:`Python` section of the user guide for more information.

**New Features and Enhancements**

Added functions to CVODE(S) and IDA(S) to set the maximum number of inequality
constraint failures in a step attempt (:c:func:`CVodeSetMaxNumConstraintFails`
and :c:func:`IDASetMaxNumConstraintFails`) and to retrieve the total number of
failed step attempts due to an inequality constraint violation
(:c:func:`CVodeGetNumConstraintFails` and
:c:func:`IDAGetNumConstraintFails`). As a result, constraint failures are no
longer included in the number of step failures due to a solver failure (i.e.,
the values returned by :c:func:`CVodeGetNumStepSolveFails` and
:c:func:`IDAGetNumStepSolveFails`). The functions
:c:func:`CVodeGetNumConstraintCorrections` and
:c:func:`IDAGetNumConstraintCorrections` were also added to retrieve the number
of steps where the corrector was modified to satisfy an inequality constraint
without failing the step.

The functions ``CVodeGetUserDataB`` and ``IDAGetUserDataB`` were added to CVODES
and IDAS, respectively.

**Bug Fixes**

Fixed a bug in the CVODE(S) inequality constraint handling where the predicted
state was used to compute the step size reduction factor which could lead to an
insufficient reduction in the step size or, when the prediction violates the
constraints, an infinitely large step size in the next step attempt (`Issue #702
<https://github.com/LLNL/sundials/issues/702>`__).

On the initial time step with a user-supplied initial step size, ARKODE and
CVODE(S) will now return ``ARK_TOO_CLOSE`` or ``CV_TOO_CLOSE``, respectively,
when the requested output time is the same as, or within numerical roundoff of,
the initial time (`Issue #722
<https://github.com/LLNL/sundials/issues/722>`__). Before a ``TOO_CLOSE`` error
would only be returned when internally estimating the initial step size. In
IDA(S), added a ``IDA_TOO_CLOSE`` return value for when the initial and output
time are too close. Previously, IDA(S) would return ``IDA_ILL_INPUT``.

Fixed a bug in ARKODE, CVODE(S), and IDA(S) where the linear solver counters
were not reset on reinitialization until the next call to advance the system. As
such, non-zero linear solver statistics could be returned if retrieving or
printing linear solver counters between reinitialization and the next call to
advance the system.

In CVODES and IDA, added missing return flag names to
:c:func:`CVodeGetReturnFlagName` and :c:func:`IDAGetReturnFlagName`,
respectively.

The SPRKStep module now accounts for zero coefficients in the SPRK tables,
eliminating extraneous function evaluations.

A bug was fixed in KINSOL where the information logging function would always be
called even when informational logging was disabled (`Issue #801
<https://github.com/LLNL/sundials/issues/801>`__).

A bug preventing a user supplied :c:func:`SUNStepper_ResetCheckpointIndex`
function from being called was fixed.

The interface to Ginkgo batched linear solvers has been updated to fix build
errors when using 64-bit index types (`Issue #797
<https://github.com/LLNL/sundials/issues/797>`__). Note, only the batched dense
matrix in Ginkgo is currently compatible with 64-bit indexing (as of Ginkgo
1.10).

The Kokkos N_Vector now properly handles unmanaged views. Previously, if a
Kokkos ``N_Vector`` was created from an unmanaged view, the view would become a
managed view and the data would be freed unexpectedly.

Fixed a CMake bug which resulted in static targets depending on shared targets
when building both types of libraries in the same build (`Issue #692
<https://github.com/LLNL/sundials/issues/692>`__).

Some installed Fortran example makefiles were not linking to
``sundials_fcore_mod`` and ``sundials_core`` libraries as they should be. This
is now fixed.

**Deprecation Notices**

The ``N_Vector_S`` typedef to ``N_Vector*`` is deprecated and will be removed in
the next major release.

The ``CSC_MAT`` and ``CSR_MAT`` macros defined in ``sunmatrix_sparse.h`` will be
removed in the next major release. Use ``SUN_CSC_MAT`` and ``SUN_CSR_MAT``
instead.

``SUNDIALSFileOpen`` and ``SUNDIALSFileClose`` will be removed in the next major
release.  Use :c:func:`SUNFileOpen` and :c:func:`SUNFileClose` instead.

The ``Convert`` methods on the ``sundials::kokkos:Vector``,
``sundials::kokkos::DenseMatrix``, ``sundials::ginkgo::Matrix``,
``sundials::ginkgo::BatchMatrix``, ``sundials::kokkos::DenseLinearSolver``,
``sundials::ginkgo::LinearSolver``, and ``sundials::ginkgo::BatchLinearSolver``
classes have been deprecated and will be removed in the next major release. The
method ``get``, should be used instead.
