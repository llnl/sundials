.. For package-specific references use :ref: rather than :numref: so intersphinx
   links to the appropriate place on read the docs

**Major Features**

**New Features and Enhancements**

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

**Deprecation Notices**
