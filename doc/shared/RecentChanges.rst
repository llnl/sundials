.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

We added a new SUNNonlinearSolver implementation,
:ref:`SUNNonlinearSolver_Auto <SUNNonlinSol.Auto>`, which uses an algorithm described in
:cite:p:`norsett1986switching` to switch between modified Newton
iteration and fixed-point iteration based on an estimate of stiffness. This
solver may be useful to pair with the BDF method in CVODE/CVODES for users who
are unsure about the stiffness of their problem.

**Bug Fixes**

Fixed a CMake bug where the SuperLU_MT interface would not be built and
installed without setting the ``SUPERLUMT_WORKS`` option to ``TRUE``.

Fixed the embedded coefficients for the ``ARKODE_TSITOURAS_7_4_5`` Butcher
table.

**Deprecation Notices**
