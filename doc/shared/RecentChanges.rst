.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Updated ``examples/cvode/petsc/cv_petsc_ex7.c`` to support PETSc 3.25.0.

Added the function :c:func:`ARKodeSkipAdaptStopTime` to specify that
stop-time-limited steps should be disregarded when selecting step sizes
for time step adaptivity.  Added the functions
:c:func:`CVodeSkipAdaptStopTime` and :c:func:`IDASkipAdaptStopTime`
to specify that stop-time-limited steps should be disregarded when
adapting the step size and method order for CVODE(S) and IDA(S), respectively.

An optional N_Vector routine, :c:func:`N_VCopy`, was added, to streamline data copies between two
vectors.  For user-supplied N_Vector modules that do not provide this function, :c:func:`N_VScale`
will be used instead.

**Bug Fixes**

Fixed memory leaks in CVODES, IDAS, and KINSOL in the unlikely event of a failed
``malloc``.

**Deprecation Notices**
