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

**Deprecation Notices**
