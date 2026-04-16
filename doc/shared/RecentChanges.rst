.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

We added a new SUNNonlinearSolver implementation,
:ref:`SUNNonlinearSolver_Auto <SUNNonlinSol.Auto>`, which uses an algorithm described in
:cite:p:`norsett1986switching` to switch between a modified Newton iteration and fixed-point
iteration based on an estimate of stiffness. This solver may be useful to pair with the BDF method
in CVODE/CVODES, or with DIRK methods in ARKODE, for users who are unsure about
the stiffness of their problem. See the module documentation for more information.

Updated ``examples/cvode/petsc/cv_petsc_ex7.c`` to support PETSc 3.25.0.

**Bug Fixes**

Fixed memory leaks in CVODES, IDAS, and KINSOL in the unlikely event of a failed
``malloc``.

**Deprecation Notices**
