..
   Author(s): David J. Gardner, Cody J. Balos @ LLNL
   -----------------------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   -----------------------------------------------------------------------------

.. _SourceCode.Naming:

Naming Conventions
==================

All exported symbols that will be publicly available must be namespaced
appropriately!

* ``SUN_`` or ``SUNDIALS_`` for macros

* ``sun`` for typedef's to native types (e.g., ``sunindextype``)

* ``SUN`` for public functions that are not class functions (see
  below for class/struct naming conventions)

* ``sun`` for private functions that are non-native.

* ``sundials::`` for C++ code (nesting under `sundials::` is OK)

* ``sundials::<somename>::impl`` for C++ code that is private (implementation
  only)

Generally Pascal case (e.g. ``DoSomething``) is used for public names and
camelcase for private names (e.g. ``doSomething``).

Macros/Constants
----------------

Upper case should be used for macros and constants.

Variable names
--------------

Snake case is preferred for local variable names e.g. ``foo_bar``.

Variables which are pointers to an array, and are effectively treated/indexed
as a contiguous array, should use the suffix `_<1|2|3>d`, e.g.

.. code-block:: c

   sunrealtype my_array[3] = {1.0, 2.0, 3.0};
   sunrealtype* sequence_1d = my_array;

   sunrealtype my_matrix[2][3] = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
   sunrealtype* my_matrix_rows[2] = { my_matrix[0], my_matrix[1] };
   sunrealtype** sequence_2d = my_matrix_rows;


Variables which are purely pointers should use the suffix, ``_ptr``, e.g.

.. code-block:: c

   N_Vector y = N_VNew_Serial(2, sunctx);
   N_Vector* y_ptr = &y;

When combining the two rules, the ``_ptr`` suffix should come last, e.g.

.. code-block:: c

   sunrealtype my_array[3] = {1.0, 2.0, 3.0};
   sunrealtype* sequence_1d = my_array;
   sunrealtype** sequence_1d_ptr = &sequence_1d;

   sunrealtype my_matrix[2][3] = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
   sunrealtype* my_matrix_rows[2] = { my_matrix[0], my_matrix[1] };
   sunrealtype** sequence_2d = my_matrix_rows;
   sunrealtype*** sequence_2d_ptr = &my_matrix_rows;


.. warning::

   The suffixes are **required** for parameters of functions within public header
   files because the Python interface generator relies on the suffixes to determine
   the proper way to expose the parameter to Python users. It is preferable to follow
   this convention within other code, but not required.


C function names
----------------

Functions should have a descriptive name that lets a reader know what it does.
For functions that return a boolean value, prefer the convention
``Is<statement>``, e.g. ``IsOutputRank``. Pascal case (with the appropriate
namespace prefix) should be used for all public function names that are not
class functions.

C++ function names
------------------

All functions should be in the ``sundials::`` namespace. Pascal case should be
used for public function names. Camelcase should be used for private function
names.

Names for Vectors, Matrices, and Solvers
----------------------------------------

The SUNDIALS vector, matrix, linear solver, and nonlinear solver classes use the
naming convention ``<short class name><method>`` for base class methods where
each component of the name uses Pascal case. See
:numref:`SourceCode.Naming.Table.OldBaseClassMethodNaming` for examples.

.. note::

   This naming convention *only* applies to the vector, matrix, and solver
   classes. All other classes should follow the naming convention described in
   :ref:`SourceCode.Naming.NewClasses`.

.. _SourceCode.Naming.Table.OldBaseClassMethodNaming:

.. Table:: SUNDIALS base class naming convention examples for vectors, matrices,
           linear solvers and nonlinear solvers.

   +-----------------------+------------------+------------+-----------------------+
   | Base Class            | Short Name       | Operation  | Method                |
   +-----------------------+------------------+------------+-----------------------+
   | ``N_Vector``          | ``N_V``          | Linear Sum | ``N_VLinearSum``      |
   +-----------------------+------------------+------------+-----------------------+
   | ``SUNMatrix``         | ``SUNMat``       | Zero       | ``SUNMatZero``        |
   +-----------------------+------------------+------------+-----------------------+
   | ``SUNLinearSolver``   | ``SUNLinSol``    | Setup      | ``SUNLinSolSetup``    |
   +-----------------------+------------------+------------+-----------------------+
   | ``SUNNonlinarSolver`` | ``SUNNonlinSol`` | Solve      | ``SUNNonlinSolSolve`` |
   +-----------------------+------------------+------------+-----------------------+

Derived class implementations of the base class methods should follow the naming
convention ``<short class name><method>_<implementation>``. See
:numref:`SourceCode.Naming.Table.OldDerivedClassMethodNaming` for examples.

.. _SourceCode.Naming.Table.OldDerivedClassMethodNaming:

.. Table:: SUNDIALS derived class naming convention examples for vectors,
           matrices, linear solvers and nonlinear solvers.

   +---------------+-----------------------+------------------------------+
   | Derived Class | Base Class Method     | Method Implementation        |
   +---------------+-----------------------+------------------------------+
   | Serial        | ``N_VLinearSum``      | ``N_VLinearSum_Serial``      |
   +---------------+-----------------------+------------------------------+
   | Dense         | ``SUNMatZero``        | ``SUNMatZero_Dense``         |
   +---------------+-----------------------+------------------------------+
   | SPGMR         | ``SUNLinSolSetup``    | ``SUNLinSolSetup_SPGMR``     |
   +---------------+-----------------------+------------------------------+
   | Newton        | ``SUNNonlinSolSolve`` | ``SUNNonlinSolSolve_Newton`` |
   +---------------+-----------------------+------------------------------+

Implementation specific methods do not currently have a consistent naming
convention across the different derived classes. When adding new methods to an
existing class, follow the naming style used within that class. When adding a
new derived class, use the same style as above for implementations of the base
class method i.e., ``<short class name><method>_<implementation>``.

.. _SourceCode.Naming.NewClasses:

Names for New Classes
---------------------

All new base classes should use the naming convention ``<class name>_<method>``
for the base class methods. See
:numref:`SourceCode.Naming.Table.NewBaseClassMethodNaming` for examples.

.. _SourceCode.Naming.Table.NewBaseClassMethodNaming:

.. Table:: SUNDIALS naming conventions for methods in new base classes.

   +-----------------------+------------+---------------------------+
   | Base Class            | Operation  | Method                    |
   +-----------------------+------------+---------------------------+
   | ``SUNMemoryHelper``   | Alloc      | ``SUNMemoryHelper_Alloc`` |
   +-----------------------+------------+---------------------------+

Derived class implementations of the base class methods should follow the naming
convention  ``<class name>_<method>_<implementation>``. See
:numref:`SourceCode.Naming.Table.NewDerivedClassMethodNaming` for examples.

.. _SourceCode.Naming.Table.NewDerivedClassMethodNaming:

.. Table:: SUNDIALS naming conventions for derived class implementations of
           methods in new base classes.

   +---------------+---------------------------+--------------------------------+
   | Derived Class | Base Class Method         | Method Implementation          |
   +---------------+---------------------------+--------------------------------+
   | CUDA          | ``SUNMemoryHelper_Alloc`` | ``SUNMemoryHelper_Alloc_Cuda`` |
   +---------------+---------------------------+--------------------------------+

For destructor functions, use ``Destroy`` rather than ``Free`` or some other alternative.


.. _SourceCode.Naming.CppClasses:

Naming Convention for C++ Classes
---------------------------------

C++ classes should have a descriptive name. The class name should not be
prefixed with ``SUN``, but it should reside in the ``sundials::`` namespace.
Public C++ class functions should use Pascal case (e.g. ``DoSomething``).
Private C++ class functions should use camelcase (e.g. ``doSomething``).

C++ private class members should use snake case with a trailing underscore
(e.g. ``some_var_``).


.. _SourceCode.Naming.Enums:

Enums
-----

Enum tags/identifiers should follow class naming rules and use Pascal case.
Enum values should follow the rules for macros and constants.
