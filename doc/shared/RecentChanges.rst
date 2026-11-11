.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added initial HIP support to sundials4py. Builds with the HIP N_Vector enabled
now expose HIP vector constructors, device pointer access, explicit host/device
copies, and the HIP memory helper.
ROCm CuPy, PyTorch, and JAX array accessors support HIP heat equation examples,
including immutable JAX arrays and JAX Refs.
Added a ``rocm`` sundials4py extra that selects a separate
``sundials4py-rocm`` wheel and requires a ROCm build of PyTorch. The ROCm wheel
loads PyTorch's HIP runtime before its native extension and rejects a
different loaded runtime.

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

Added support for the MPIPlusX and MPIManyVector ``N_Vector`` modules to
sundials4py in source builds with ``-DSUNDIALS_ENABLE_MPI=ON``. Communicators
are passed as ``mpi4py`` communicators, and ``None`` (and, in all builds, any
integer-like object equal to zero other than a ``bool``) is now accepted as a
synonym for ``SUN_COMM_NULL``.

**Bug Fixes**

Fixed a bug in the sundials4py CUDA ``N_Vector`` DLPack interface where the
DLPack capsule always reported device 0 (``__dlpack_device__`` was correct), so
PyTorch, CuPy and JAX arrays of vectors on other CUDA devices were mislabeled.

Fixed a segfault that could occur when handling errors without a ``SUNContext``
provided.

Fixed a bug in ``FindMAGMA.cmake`` which didn't allow use of MAGMA versions with
multiple digits in an identifier.

Fixed a bug in the LSRKStep module where the step post-processing function for
the RKL step was called with an incorrect time argument.

Fixed a bug in :c:func:`LSRKStepSetDomEigFrequency` where a negative ``nsteps``
value was stored as the update frequency instead of restoring the default value.

**Deprecation Notices**
