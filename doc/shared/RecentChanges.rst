.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

SUNDIALS now has official Python interfaces! With this release, we are shipping a **beta version** of
the sundials4py Python module (created with nanobind and litgen). sundials4py provides explicit
interfaces to most features of SUNDIALS. 

**New Features and Enhancements**

The functions ``CVodeGetUserDataB`` and ``IDAGetUserDataB`` were added to CVODES
and IDAS, respectively.

**Bug Fixes**

Fixed a CMake bug which resulted in static targets depending on shared targets
when building both types of libraries in the same build
(`Issue #692 <https://github.com/LLNL/sundials/issues/692>`__).

Some installed Fortran example makefiles were not linking to ``sundials_fcore_mod``
and ``sundials_core`` libraries as they should be. This is now fixed.

On the initial time step with a user-supplied initial step size, ARKODE and
CVODE(S) will now return ``ARK_TOO_CLOSE`` or ``CV_TOO_CLOSE``, respectively,
when the requested output time is the same as the initial time (or within
numerical roundoff of the initial time). Before a ``TOO_CLOSE`` error would only
be returned when internally estimating the initial step size. In IDA(S), added a
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

A bug was fixed in KINSOL where the information logging function would always be
called even when informational logging was disabled.

**Deprecation Notices**

``SUNDIALSFileOpen`` and ``SUNDIALSFileClose`` will be removed in the next major release. 
Use :c:func:`SUNFileOpen` and :c:func:`SUNFileClose` instead.

The ``Convert`` methods on the ``sundials::kokkos:Vector``, ``sundials::kokkos::DenseMatrix``,
``sundials::ginkgo::Matrix``, ``sundials::ginkgo::BatchMatrix``, ``sundials::kokkos::DenseLinearSolver``,
``sundials::ginkgo::LinearSolver``, and ``sundials::ginkgo::BatchLinearSolver`` classes have
been deprecated and will be removed in the next major release. The method ``get``, should
be used instead.

The ``CSC_MAT`` and ``CSR_MAT`` macros defined in ``sunmatrix_sparse.h`` will be removed in
the next major release. Use ``SUN_CSC_MAT`` and ``SUN_CSR_MAT`` instead.

The ``N_Vector_S`` typedef to ``N_Vector*`` is deprecated and will be removed in the next major release.
