.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added the function :c:func:`SUNLogger_SetQueueAndFlushMsgFns` to allow for
user-defined functions to queue and flush log messages.

Updated ``examples/cvode/petsc/cv_petsc_ex7.c`` to support PETSc 3.25.0.

Added the ``ARKODE_SSP_ERK_3_1_2``, ``ARKODE_SSP_ERK_4_1_2``,
``ARKODE_SSP_ERK_4_2_3``, ``ARKODE_SSP_ERK_10_3_4``, ``ARKODE_SSP_LSPUM_ERK_3_1_2``,
and ``ARKODE_ASCHER_ERK_3_1_2``
embedded explicit Runge-Kutta Butcher tables.

Added the ``ARKODE_SSP_DIRK_3_1_2``, ``ARKODE_SSP_LSPUM_SDIRK_3_1_2``,
``ARKODE_SSP_ESDIRK_4_2_3``, and ``ARKODE_ASCHER_SDIRK_3_1_2`` embedded diagonally
implicit Runge-Kutta Butcher tables.

Of these, embedded additive Runge-Kutta methods may be formed using
``ARKODE_SSP_ERK_3_1_2`` + ``ARKODE_SSP_DIRK_3_1_2``,
``ARKODE_SSP_ERK_4_2_3`` + ``ARKODE_SSP_ESDIRK_4_2_3``,
``ARKODE_SSP_LSPUM_ERK_3_1_2`` + ``ARKODE_SSP_LSPUM_SDIRK_3_1_2``,
and ``ARKODE_ASCHER_ERK_3_1_2`` + ``ARKODE_ASCHER_SDIRK_3_1_2``.

Added the ``ARKODE_IMEX_MRI_GARK_ARS222`` and ``ARKODE_IMEX_MRI_GARK_GIRALDO``
embedded implicit-explicit MRI-GARK coupling tables.

**Bug Fixes**

Fixed memory leaks in CVODES, IDAS, and KINSOL in the unlikely event of a failed
``malloc``.

Fixed minor bug in reporting the maximum number of stages in
:c:func:`ARKodeGetStageIndex` when running SSP methods in LSRKStep.

Removed duplicate logging output that would cause the Python logging tools to
fail with a repeated key error.

Fixed a CMake issue that prevented finding third-party libraries installed in
default search locations e.g., paths included in ``CMAKE_INSTALL_PREFIX``
(`Issue #935 <https://github.com/llnl/sundials/issues/935>`__).

Fixed empty ``elseif()`` cases in the CMake files for the Fortran interfaces to
the ManyVector and MPIPlusX vectors which could results in a missing include
path when compiling if an MPI compiler wrapper is not found.

Fixed a bug where IDAS would incorrectly compute the quadrature predictor when
IDACalcIC was used. In some cases, this lead to an inconsistent solution in the
forward solve compared to the forward recomputation from a checkpoint,
ultimately causing a segfault.

**Deprecation Notices**
