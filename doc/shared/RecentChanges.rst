.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added the utility function, :c:func:`SUNFileFlush` for flushing file
pointers. This is useful when using the Fortran 2003 interfaces.

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

Fixed a segfault that could occur when handling errors without a ``SUNContext``
provided.

Fixed a bug in ``FindMAGMA.cmake`` which didn't allow use of MAGMA versions with 
multiple digits in an identifier.

Fixed a bug in the LSRKStep module where the step post-processing function for
the RKL step was called with an incorrect time argument.

**Deprecation Notices**
