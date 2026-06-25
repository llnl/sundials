.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**New Features and Enhancements**

We added a new SUNNonlinearSolver implementation,
:ref:`SUNNonlinearSolver_Auto <SUNNonlinSol.Auto>`, which uses an algorithm described in
:cite:p:`norsett1986switching` to switch between a modified Newton iteration and fixed-point
iteration based on an estimate of stiffness. This solver may be useful to pair with the BDF method
in CVODE/CVODES, or with DIRK methods in ARKODE, for users who are unsure about
the stiffness of their problem. See the module documentation for more information. We also 
extended the :ref:`SUNNonlinearSolver API <SUNNonlinSol.API>` with callback setters
:c:func:`SUNNonlinSolSetNormFn`, :c:func:`SUNNonlinSolSetGetUpdateNormFn`, and
:c:func:`SUNNonlinSolSetGetConvRateFn`.

Added the function :c:func:`SUNLogger_SetQueueAndFlushMsgFns` to allow for
user-defined functions to queue and flush log messages.

Updated ``examples/cvode/petsc/cv_petsc_ex7.c`` to support PETSc 3.25.0.

Added the ``ARKODE_SSP_ERK_3_1_2``, ``ARKODE_SSP_ERK_4_1_2``,
``ARKODE_SSP_ERK_4_2_3``, ``ARKODE_SSP_ERK_10_3_4``, ``ARKODE_SSP_LSPUM_ERK_3_1_2``,
and ``ARKODE_ASCHER_ERK_3_1_2``
embedded explicit Runge-Kutta Butcher tables.

Added the ``ARKODE_SSP_DIRK_3_1_2``, ``ARKODE_SSP_LSPUM_SDIRK_3_1_2``,
``ARKODE_ESDIRK_4_2_3``, and ``ARKODE_ASCHER_SDIRK_3_1_2`` embedded diagonally
implicit Runge-Kutta Butcher tables.

Of these, embedded additive Runge-Kutta methods may be formed using
``ARKODE_SSP_ERK_3_1_2`` + ``ARKODE_SSP_DIRK_3_1_2``,
``ARKODE_SSP_ERK_4_2_3`` + ``ARKODE_ESDIRK_4_2_3``,
``ARKODE_SSP_LSPUM_ERK_3_1_2`` + ``ARKODE_SSP_LSPUM_SDIRK_3_1_2``,
and ``ARKODE_ASCHER_ERK_3_1_2`` + ``ARKODE_ASCHER_SDIRK_3_1_2``.

Added the ``ARKODE_IMEX_MRI_GARK_ASCHER_ARK2`` and ``ARKODE_IMEX_MRI_GARK_ARK2``
embedded implicit-explicit MRI-GARK coupling tables.

When info or debug logging is enabled, i.e., ``SUNDIALS_LOGGING_LEVEL`` is at
least 3, it will now print to ``stdout`` by default. Previously, the default was
to produce no output, and this behavior can be restored by setting the
environment variables ``SUNLOGGER_INFO_FILENAME`` and
``SUNLOGGER_DEBUG_FILENAME`` to an empty string.

Improved the performance of logging when enabled but no file pointer was set.

**Bug Fixes**

Fixed a bug in the sundials4py wrappers for user-provided SUNStepper evolve
and one_step functions.

Fixed a bug where an unrecognized error return flag from a user-provided
SUNLinearSolver module would register as a successful linear solve.

Fixed a minor bug where the number of required stages for STS methods
in the LSRKStep module was incorrectly computed using the spectral
radius instead of the real part of the Jacobian eigenvalues.

Fixed a minor bug where the negative real extent of the stability region
for the RKC method was not being properly computed, which could result
in an underestimation of the number of stages required for stability.

Fixed a minor bug where STS methods were limited to one fewer than
the maximum allowed number of stages. STS can now use the full maximum
number of stages.

Fixed a bug that caused SUNDomEigEstimator_Initialize to overwrite the
user-provided initial guess from SUNDomEigEstimator_SetInitialGuess. These
routines are now order-independent.

Fixed memory leaks in CVODES, IDAS, and KINSOL in the unlikely event of a failed
``malloc``.

Fixed minor bug in reporting the maximum number of stages in
:c:func:`ARKodeGetStageIndex` when running SSP methods in LSRKStep.

Removed duplicate logging output that would cause the Python logging tools to
fail with a repeated key error.

Fixed a CMake issue that prevented finding third-party libraries installed in
default search locations e.g., paths included in ``CMAKE_INSTALL_PREFIX``
(`Issue #935 <https://github.com/llnl/sundials/issues/935>`__).

Fixed a CMake issue that prevented automatically finding PETSc dependencies.

Fixed empty ``elseif()`` cases in the CMake files for the Fortran interfaces to
the ManyVector and MPIPlusX vectors which could results in a missing include
path when compiling if an MPI compiler wrapper is not found.

Fixed a bug where IDAS would incorrectly compute the quadrature predictor when
IDACalcIC was used. In some cases, this lead to an inconsistent solution in the
forward solve compared to the forward recomputation from a checkpoint,
ultimately causing a segfault.

Fixed a CMake bug where Fortran modules were not created for LSRKStep,
ForcingStep, and SplittingStep.

Corrected the version number used in version added, changed, and deprecated
notes in the documentation to always use the SUNDIALS version number with the
package version number as a parenthetical note when it differs from the SUNDIALS
version number.

Fixed a bug in sundials4py where the :c:func:`CVodeGetRootInfo`,
:c:func:`ARKodeGetRootInfo`, and :c:func:`IDAGetRootInfo` functions did not
correctly return the rootsfound array. This addresses `Issue #937
<https://github.com/llnl/sundials/issues/937>`__.

Fixed a bug in ERKStep where calling :c:func:`ARKodeResize` before
:c:func:`ARKodeEvolve` or :c:func:`ARKodeInit` would result in a segmentation
fault.
