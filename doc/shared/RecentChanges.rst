.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

Added the utility function, :c:func:`SUNFileFlush` for flushing file
pointers. This is useful when using the Fortran 2003 interfaces.

sundials4py now supports implementing SUNDIALS classes in Python. The new
``CustomSUNMatrix``, ``CustomSUNLinearSolver``, ``CustomSUNNonlinearSolver``,
``CustomSUNHController``, and ``CustomSUNMRIController`` base classes may be
subclassed to provide a :c:type:`SUNMatrix`, :c:type:`SUNLinearSolver`,
:c:type:`SUNNonlinearSolver`, or :c:type:`SUNAdaptController` implementation
written in Python, and instances of such a subclass may be passed to any
sundials4py function that expects the corresponding SUNDIALS object. See
:ref:`Python.Usage.CustomObjects` for details, and the
``cvs_custom_nonlinsol.py``, ``kin_custom_linsol.py``, and
``ark_custom_adaptcontroller.py`` examples for annotated templates.

**Bug Fixes**

Fixed a segfault that could occur when handling errors without a ``SUNContext``
provided.

Fixed a bug in ``FindMAGMA.cmake`` which didn't allow use of MAGMA versions with 
multiple digits in an identifier.

**Deprecation Notices**
