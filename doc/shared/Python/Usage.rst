.. ----------------------------------------------------------------
   SUNDIALS Copyright Start
   Copyright (c) 2025-2026, Lawrence Livermore National Security,
   University of Maryland Baltimore County, and the SUNDIALS contributors.
   Copyright (c) 2013-2025, Lawrence Livermore National Security
   and Southern Methodist University.
   Copyright (c) 2002-2013, Lawrence Livermore National Security.
   All rights reserved.

   See the top-level LICENSE and NOTICE files for details.

   SPDX-License-Identifier: BSD-3-Clause
   SUNDIALS Copyright End
   ----------------------------------------------------------------

.. _Python.Usage:

Using sundials4py
=================

At a high level, using SUNDIALS from Python via sundials4py looks a lot like
using SUNDIALS from C or C++. Below we overview using sundials4py and discuss
the few notable differences.

.. _Python.Usage.Installation:

Installation
------------

You can install sundials4py directly from `PyPI
<https://pypi.org/project/sundials4py/>`__ using pip:

.. code-block:: bash

   pip install sundials4py

You can also install sundials4py from git:

.. code-block:: bash

   pip install git+https://github.com/LLNL/sundials.git

The default build of sundials4py that is distributed as a binary wheel uses
double precision real types and 64-bit indices. To install SUNDIALS with
different precisions and index sizes, you can build from source wheels instead
of using the pre-built binary wheels. When building from source wheels instead
of binary wheels, you can customize the SUNDIALS precision (real type) and index
type at build time by passing the CMake arguments in an environment variable
when running pip. For example:

.. code-block:: bash

   export CMAKE_ARGS="-DSUNDIALS_PRECISION=SINGLE -DSUNDIALS_INDEX_SIZE=64"
   pip install sundials4py --no-binary=sundials4py

Other SUNDIALS options can also be accessed in this way. Review
:numref:`Installation.Options` for more information on the available options.

CUDA support
^^^^^^^^^^^^

CUDA support is optional and must be enabled when building sundials4py from
source. The CUDA N_Vector module is enabled by default when CUDA support is
enabled. For example:

.. code-block:: bash

   export CMAKE_ARGS="-DSUNDIALS_ENABLE_CUDA=ON"
   pip install sundials4py --no-binary=sundials4py

With a CUDA-enabled build, the Python interface supports the following array
backends for CUDA N_Vectors:

- CuPy device views through ``N_VGetCupyArray``
- PyTorch native-device tensor views through ``N_VGetTorchTensor``
- JAX native-device array views through ``N_VGetJaxArray``

The required array framework must be installed separately. CPU-only builds
continue to support NumPy and the CPU modes of PyTorch and JAX.

KLU support
^^^^^^^^^^^

To build sundials4py with KLU support enabled, pass CMake options for KLU (see :numref:`Installation.Options.KLU`) through ``CMAKE_ARGS`` when building from
source:

.. code-block:: bash

   export CMAKE_ARGS="-DSUNDIALS_ENABLE_KLU=ON -DKLU_ROOT=/path/to/suitesparse/installation"
   pip install sundials4py --no-binary=sundials4py

.. _Python.Usage.Modules:

Modules
-------

After installation, you can import the sundials4py module with

.. code-block:: python

   import sundials4py

which includes the following submodules (which may also be individually
imported) for accessing specific SUNDIALS features:

- ``sundials4py.core`` contains all the shared SUNDIALS classes and functions as
  well as many of the native SUNDIALS class implementations:

  - NVector: serial and many-vector

  - SUNMatix: band, dense, and sparse

  - SUNLinearSover: band, dense, PCG, SPBCGS, SPFGMR, SPGMR, and SPTFQMR

  - SUNNonlinearSolver: fixed-point and Newton

  - SUNAdaptController: Soderlind, ImEx-Gus, and MRI H-Tol

  - SUNDomEigEstimator: Power

  - SUNAdjointCheckPointScheme: Fixed

- ``sundials4py.arkode`` contains all of the ARKODE specific classes and
  functions

- ``sundials4py.cvodes`` contains all of the CVODES specific classes and
  functions

- ``sundials4py.idas`` contains all of the IDAS specific classes and functions

- ``sundials4py.kinsol`` contains all of the KINSOL specific classes and
  functions

CVODE and IDA dot not have modules because CVODES and IDAS provide all of the
same capabilities plus continuous forward and adjoint sensitivity analysis.

.. note::

   Not all SUNDIALS features are supported by the Python interfaces. In
   particular, third-party libraries are only partially supported.

.. _Python.Usage.Example:

Example Usage
-------------

We now consider a simple CVODE example to illustrate using sundials4py and
highlight some of the differences to using SUNDIALS from C/C++. The items
highlighted below similarly apply to using other SUNDIALS packages. For more
information on usage differences, continue to the :ref:`next section
<Python.Usage.Differences>`. Additional examples can be found in the
``examples/python`` directory of the :examples:`SUNDIALS GitHub repository <python>`.

This example demonstrates how to use CVODES to solve the Lotka-Volterra
equations, a model of predator-prey dynamics in ecology, given by

.. math::

   u' &=  p_0 u - p_1 u v \\
   v' &= -p_2 v + p_3 u v

where :math:`u` is the prey population, :math:`v` is the predator population,
:math:`p_0` is prey birth rate, :math:`p_1` is the predation rate, :math:`p_2`
is the predator death rate, and :math:`p_3` is predator growth rate from
predation. We use the parameters :math:`p = [1.5, 1.0, 3.0, 1.0]`, initial
condition :math:`y(0) = [1.0, 1.0]`, and integration interval :math:`t \in [0,
10]`.

.. literalinclude:: cvs_lotkavolterra.py
   :language: python
   :start-after: # --- start example ---
   :end-before: # --- end example ---
   :linenos:
   :emphasize-lines: 4-5,36,49,74,108,119-122,149,250-254

.. _Python.Usage.Differences:

Usage Differences
-----------------

While sundials4py closely follows the C API, some differences are inevitable due
to the differences between Python and C as well as the requirements of the code
generation tool used to create the bindings. In this section, we note the most
critical differences.


View Classes and Memory Management
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

sundials4py provides natural usage of SUNDIALS objects with object lifetimes
managed by the Python garbage collection as with any other Python object. There
is only one caveat, the SUNDIALS integrator/solver ``void*`` objects (those
returned by ARKODE, CVODES, IDAS, and KINSOL ``Create`` constructors) are
wrapped in "View" classes (behind the scenes) for compatibility with
nanobind. These view objects cannot be implicitly converted to the underlying
``void*``. As such, when calling a function which operates on these ``void*``
objects, one must extract the ``void*`` "capsule" from the view object by
calling the view's ``get`` method.

.. code-block:: python

   # Create CVODE object (returns void* in C)
   cvode = CVodeCreate(CV_BDF, sunctx)

   # Notice we need to call cvode.get()
   status = CVodeInit(cvode.get(), ode_problem.f, T0, y)


.. _Python.Usage.Differences.Lifetimes:

Lifetime of Objects Retained by SUNDIALS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Some SUNDIALS functions retain pointers to objects supplied by the
application, such as matrices, linear and nonlinear solvers, and adaptivity
controllers. The Python application must keep each supplied object alive for as
long as the SUNDIALS integrator or solver may use it. This requirement applies
to both native sundials4py objects and Python subclasses of the custom SUNDIALS
object interfaces.

Store retained objects in variables whose lifetime covers all SUNDIALS calls
that may use them. Do not rely on a temporary object passed directly to a
setter, and do not delete or overwrite the final Python reference before the
receiving SUNDIALS object is destroyed or reconfigured.

.. code-block:: python

   A = MyCustomMatrix(..., sunctx)
   LS = MyCustomLinearSolver(..., sunctx)

   # CVODE retains the native pointers associated with LS and A.
   CVodeSetLinearSolver(cvode_mem.get(), LS, A)

   # Keep LS and A alive while cvode_mem may use them.
   CVode(cvode_mem.get(), tout, y, tret, CV_NORMAL)


Return-by-Pointer Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Functions that return values via pointer arguments in the C API are mapped to
Python functions that return a tuple where the **first element** is the
function's return value (typically an error code) and **subsequent elements**
are the values that would be returned via pointer arguments in C, in the same
order as the C function signature.

**Example 1: Single Return-by-Pointer Value**

C:
   .. code-block:: C

      int retval;
      long int numsteps;
      retval = CVodeGetNumSteps(cvode_mem, &numsteps);
      printf("Number of steps: %ld\n", numsteps);

Python:
   .. code-block:: python

      retval, numsteps = CVodeGetNumSteps(cvode_mem.get())
      print(f"Number of steps: {numsteps}")

**Example 2: Multiple Return-by-Pointer Values**

C:
   .. code-block:: C

      int retval;
      long int nni, ncfn;
      retval = CVodeGetNonlinSolvStats(cvode_mem, &nni, &ncfn)
      printf("Nonlinear iterations: %ld, Nonlinear convergence fails: %ld\n", nni, ncfn);

Python:
   .. code-block:: python

      retval, nni, ncfn = CVodeGetNonlinSolvStats(cvode_mem.get())
      print(f"Nonlinear iterations: {nni}, Nonlinear convergence fails: {ncnf}");


Arrays
^^^^^^

CPU ``N_Vector`` objects in sundials4py are compatible with numpy's
``ndarray``. Each CPU ``N_Vector`` can work on a numpy array without copies,
and you can access and modify the underlying data directly using
``N_VGetNumpyArray``, which returns a numpy ``ndarray`` view of the data.

SUNDIALS matrix types (dense, banded, sparse) are also exposed as Python objects
that provide access to their underlying data as numpy arrays (e.g., via
``SUNDenseMatrix_Data``).

Arrays of scalars (e.g., scaling factors passed to ``N_VLinearCombination``) are
also represented as numpy arrays.

**Example: Accessing and modifying an N_Vector**

.. code-block:: python

   y_nvec = N_VNew_Serial(10, sunctx)
   y = N_VGetNumpyArray(y_nvec)
   y[:] = np.linspace(0, 1, 10)  # Set values using numpy

**Example: Using a matrix as a numpy array**

.. code-block:: python

   mat = SUNDenseMatrix(3, 3, sunctx)
   arr = SUNDenseMatrix_Data(mat)
   arr[:] = np.eye(3)  # Set to identity matrix

This allows you to use numpy operations for vector and matrix data, and to pass
numpy arrays to and from SUNDIALS routines efficiently and without unnecessary
copies.

For CUDA ``N_Vector`` objects, ``N_VGetNumpyArray`` is not supported. Use the
backend-specific accessor instead. For example, a PyTorch CUDA tensor can be
converted to NumPy explicitly with:

.. code-block:: python

   tensor = N_VGetTorchTensor(y_nvec)
   host_array = tensor.cpu().numpy()

Similarly, use ``cupy_array.get()`` for CuPy or
``np.asarray(N_VGetJaxArray(y_nvec))`` for JAX.

For CPU ``N_Vector`` objects, ``N_VGetNumpyArray`` and
``N_VGetTorchTensor`` create lightweight zero-copy views. ``N_VGetJaxArray``
also exposes the vector data without copying, but importing an external buffer
into the JAX runtime has appreciably higher fixed overhead. In performance
sensitive CPU callbacks or Python loops, obtain a JAX view once and reuse it
rather than creating a new view for every operation. JAX is generally most
effective when its JIT compilation or GPU execution amortizes this setup cost.

The CUDA Python interface currently treats ``cuda`` as a device type rather
than accepting a CUDA device ordinal. A CUDA ``N_Vector`` uses the device
selected by its CUDA configuration. Support for multiple CUDA devices is a
future enhancement.

An ``N_Vector`` is one-dimensional by definition, so the Python accessors
require one-dimensional arrays. Applications that need a multidimensional
view can reshape the returned contiguous array in Python while keeping the
underlying vector storage one-dimensional.


User-Supplied Callback Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SUNDIALS packages and several modules/classes require user-supplied callback
functions to define problem-specific behavior, such as the right-hand side of an
ODE or a nonlinear system function. In sundials4py, you can provide these as
standard Python functions or lambdas.

The callback signatures follow the C API. As such, ``N_Vector`` arguments are
passed as ``N_Vector`` objects and the underlying ndarray must be extracted in
the user code. The only caveat is that return-by-pointer parameters are removed
from the signature, and instead become return values (mirroring how
return-by-pointer parameters for other functions are handled)

.. warning::

   The C function signatures for most callbacks include a ``void* user_data``
   argument. In Python, this argument must be present in the signature, but it
   should be ignored to avoid catastrophic errors. We recommend using ``_``
   as the parameter name in the callback signature to indicate this argument
   is unused.

**Example: ODE right-hand side for ARKStep**

.. code-block:: python

   # The C signature is:
   # int(sunrealtype t, N_Vector y, N_Vector ydot, void* user_data)
   def rhs(t, y_nvector, ydot_nvector, _): # note _ in place of user_data
      # Compute ydot = f(t, y)
      y = N_VGetNumpyArray(y_nvector)
      ydot = N_VGetNumpyArray(ydot_nvector)
      ydot[:] = -y
      return 0

   ark = ARKStepCreate(rhs, None, t0, y, sunctx)

**Example: Nonlinear system for KINSOL**

.. code-block:: python

   # The C signature is:
   # int(N_Vector u, N_Vector g, void* user_data)
   def fp_function(u_nvector, g_nvector, _): # note _ in place of user_data
      # Compute g = F(u)
      u = N_VGetNumpyArray(u_nvector)
      g = N_VGetNumpyArray(g_nvector)
      g[:] = u**2 - 1
      return 0

   kin = KINCreate(sunctx)
   KINInit(kin.get(), fp_function, u)

**Example: ARKODE LSRKStep dominant eigenvalue estimation function with return-by-pointer parameters**

.. code-block:: python

   # The C signature is:
   # int(sunrealtype t, N_Vector y, N_Vector fn,
   #     sunrealtype* lambdaR, sunrealtype* lambdaI,
   #     void* user_data,
   #     N_Vector temp1, N_Vector temp2, N_Vector temp3)
   def dom_eig(t, yvec, fnvec, _, temp1, temp2, temp3): # note the _ in place of user_data
        lamdbaR = L
        lamdbaI = 0.0
        # lambdaR and lambdaI should be returned in the order that they appear
        # as parameters in the C API and follow the error code to return
        return 0, lamdbaR, lamdbaI


Error Codes
^^^^^^^^^^^

The named ``SUN_ERR_*`` code constants are not available in Python. However, all
negative values of ``SUNErrCode`` are still errors, zero is success, and
positive values are warnings. As such, users Users can call ``SUNGetErrMsg``
from Python with the returned ``SUNErrCode`` to get further information about an
error.


.. _Python.Usage.CustomObjects:

Implementing SUNDIALS Objects in Python
---------------------------------------

SUNDIALS is designed so that applications may supply their own implementations of
several of its class interfaces. sundials4py exposes four of those interfaces for
subclassing directly in Python, so that a matrix, linear solver, nonlinear
solver, or time step controller written in Python can be handed to any SUNDIALS
package exactly as a C implementation would be:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Python base class
     - C class
     - Required method(s)
   * - ``CustomSUNMatrix``
     - :c:type:`SUNMatrix`
     - ``clone``, ``zero``, ``copy``, ``scaleadd``, ``scaleaddi``, ``matvec``
   * - ``CustomSUNLinearSolver``
     - :c:type:`SUNLinearSolver`
     - ``solve``
   * - ``CustomSUNNonlinearSolver``
     - :c:type:`SUNNonlinearSolver`
     - ``solve``
   * - ``CustomSUNHController``
     - :c:type:`SUNAdaptController` (``SUN_ADAPTCONTROLLER_H``)
     - ``estimate_step``
   * - ``CustomSUNMRIController``
     - :c:type:`SUNAdaptController` (``SUN_ADAPTCONTROLLER_MRI_H_TOL``)
     - ``estimate_step_tol``

All six classes are provided by the ``sundials4py.core`` module.
``CustomSUNHController`` and ``CustomSUNMRIController`` share the base class
``CustomSUNAdaptController``, which is not intended to be subclassed directly.

Three annotated templates are provided as a starting point:

* ``examples/python/kinsol/kin_custom_linsol.py`` -- a ``SUNMatrix`` and a
  ``SUNLinearSolver`` driving KINSOL
* ``examples/python/cvodes/cvs_custom_nonlinsol.py`` -- a
  ``SUNNonlinearSolver`` driving CVODES
* ``examples/python/arkode/ark_custom_adaptcontroller.py`` -- a
  ``SUNAdaptController`` driving ARKODE


.. _Python.Usage.CustomObjects.HowItWorks:

How Subclassing Works
^^^^^^^^^^^^^^^^^^^^^

To implement one of these interfaces, subclass the corresponding base class,
override the SUNDIALS operations you wish to support, and call the base class
constructor. The constructor takes the ``SUNContext``, plus a solver type for
the two solver families:

.. code-block:: python

   import numpy as np
   from sundials4py.core import *

   class MyMatrix(CustomSUNMatrix):
       def __init__(self, n, sunctx):
           self.data = np.zeros((n, n))
           # Call the base constructor LAST: it makes the object convertible to
           # a native handle, so the operations below must be ready to run.
           super().__init__(sunctx)

       def zero(self):
           self.data[:] = 0.0
           return SUN_SUCCESS

       # ... the remaining required operations ...

An instance is an ordinary Python object; the native SUNDIALS handle is created
lazily, the first time a binding needs one, and is then cached for the lifetime
of the Python object. Instances are passed to SUNDIALS functions directly, with
no explicit conversion step:

.. code-block:: python

   A = MyMatrix(n, sunctx)
   status = SUNMatZero(A)                      # handle created here
   status = KINSetLinearSolver(kin.get(), LS, A)

The required operations listed in the table above are validated when the native
handle is created, not at construction time, so constructing an incomplete
subclass succeeds and the first SUNDIALS function it is passed to raises
``TypeError`` instead. The message is the usual nanobind argument-conversion
error, which reports the function that rejected the object but not the reason,
so a ``TypeError`` naming a function that was given a subclass of one of these
classes almost always means a required operation is missing.

Every other operation of the interface is optional, and is detected by
inspecting the subclass rather than by any registration step: overriding a
method installs it in the native operation table, and leaving it alone leaves
that entry empty, so SUNDIALS behaves exactly as it does for a C implementation
that does not provide the operation. An intermediate base class shared by
several implementations counts as an override, so a mixin may supply operations.


.. _Python.Usage.CustomObjects.Conventions:

Method Conventions
^^^^^^^^^^^^^^^^^^

The Python methods follow the C operations closely, with the sundials4py
conventions described in :ref:`Python.Usage.Differences`:

* Methods return the SUNDIALS status code, and operations with
  return-by-pointer outputs return a tuple of the status code followed by
  those outputs, in the order they appear in the C signature. For example,
  ``estimate_step`` returns ``(status, hnew)``.

* ``N_Vector`` arguments arrive as ``N_Vector`` objects, not as arrays. Reach
  the data with ``N_VGetNumpyArray``, or express the operation with the ``N_V*``
  functions, rather than assuming a particular vector implementation.

* Methods that receive another object of the *same* class -- such as
  ``copy(dst)`` and ``scaleadd(c, other)`` on a matrix -- receive the Python
  object, so your own attributes are available on it.

.. warning::

   Operations that receive an object of a *different* class receive it as an
   opaque native handle, not as the Python object behind it. In particular, the
   ``A`` argument to a linear solver's ``setup(A)`` and ``solve(A, x, b, tol)``
   is a native :c:type:`SUNMatrix`, even when the matrix is a
   ``CustomSUNMatrix`` subclass. A Python linear solver should therefore be
   given a reference to its matrix when it is constructed, and read the entries
   through that reference. This is not usually a restriction in practice, since
   a matrix and the linear solver that factors it must already agree on how the
   entries are stored, and so are written together.

.. warning::

   Do not raise an exception to report a recoverable failure. An exception that
   escapes one of these methods is converted to the SUNDIALS error code
   ``SUN_ERR_EXT_FAIL``, which is unrecoverable, and its type, message, and
   traceback are reported through the object's :c:type:`SUNContext`. Return a
   positive status code instead: a singular matrix or a nonconverged iteration
   reported that way lets the calling package retry with a smaller step,
   whereas an exception aborts the integration.

The lifetime rules in :ref:`Python.Usage.Differences.Lifetimes` apply to these
objects: the native handle holds only a weak reference back to the Python
object, so the application must keep the Python object alive for as long as
SUNDIALS may use it. The one exception is a handle that SUNDIALS creates and
owns itself, such as the result of :c:func:`SUNMatClone`; those hold a strong
reference, so a cloned matrix keeps its Python implementation alive until
SUNDIALS destroys it.


.. _Python.Usage.CustomObjects.SUNMatrix:

CustomSUNMatrix
^^^^^^^^^^^^^^^

``CustomSUNMatrix(sunctx)`` implements the :c:type:`SUNMatrix` class. All six
operations below are required, since a calling package has no way to ask which
of them are missing. KINSOL uses only ``clone`` and ``zero``, an implicit ODE
integrator additionally needs ``scaleaddi`` to form :math:`I - \gamma J`, and an
iterative linear solver needs ``matvec``.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Method
     - Description
   * - ``clone()``
     - Return a new, empty matrix with the same shape and structure as ``self``,
       *not* a copy of the entries. Implements :c:func:`SUNMatClone`.
   * - ``zero()``
     - Set every entry of ``self`` to zero. Implements :c:func:`SUNMatZero`.
   * - ``copy(dst)``
     - Copy ``self`` into the matrix ``dst``, which is a Python object of your
       own class. Implements :c:func:`SUNMatCopy`.
   * - ``scaleadd(c, other)``
     - In place, ``self = c*self + other``, where ``other`` is a Python object
       of your own class. Implements :c:func:`SUNMatScaleAdd`.
   * - ``scaleaddi(c)``
     - In place, ``self = c*self + I``; this adds to the diagonal rather than
       overwriting it. Implements :c:func:`SUNMatScaleAddI`.
   * - ``matvec(x, y)``
     - Compute ``y = self*x`` for ``N_Vector`` arguments ``x`` and ``y``.
       Implements :c:func:`SUNMatMatvec`.

The following are optional:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Method
     - Description
   * - ``matvecsetup()``
     - Perform any setup needed before a sequence of ``matvec`` calls.
       Implements :c:func:`SUNMatMatvecSetup`.
   * - ``hermitian_transpose_matvec(x, y)``
     - Compute ``y = self^H * x``. Implements
       :c:func:`SUNMatHermitianTransposeVec`.


.. _Python.Usage.CustomObjects.SUNLinearSolver:

CustomSUNLinearSolver
^^^^^^^^^^^^^^^^^^^^^

``CustomSUNLinearSolver(sunctx, solver_type)`` implements the
:c:type:`SUNLinearSolver` class. ``solver_type`` is one of
``SUNLINEARSOLVER_DIRECT``, ``SUNLINEARSOLVER_ITERATIVE``, or
``SUNLINEARSOLVER_MATRIX_ITERATIVE``, and tells the calling package how to use
the solver: a direct solver is asked for an exact solve and its ``tol``
argument is ignored, while an iterative solver is expected to honor ``tol`` and
to report ``num_iters`` and ``res_norm``.

Only ``solve`` is required.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Method
     - Description
   * - ``solve(A, x, b, tol)``
     - **Required.** Solve ``A*x = b`` for the ``N_Vector`` ``x``. ``A`` is an
       opaque native :c:type:`SUNMatrix` handle (see the warning in
       :ref:`Python.Usage.CustomObjects.Conventions`). Implements
       :c:func:`SUNLinSolSolve`.
   * - ``initialize()``
     - Perform any setup that must happen after construction but before first
       use. Implements :c:func:`SUNLinSolInitialize`.
   * - ``setup(A)``
     - Prepare for subsequent ``solve`` calls, for example by factoring the
       matrix. Called whenever the package believes the matrix has changed;
       ``solve`` may be called many times before the next ``setup``. Implements
       :c:func:`SUNLinSolSetup`.
   * - ``set_atimes(atimes)``
     - Receive the matrix-vector product callback as a Python callable
       ``atimes(x, y)``, for a matrix-free iterative solver. Implements
       :c:func:`SUNLinSolSetATimes`.
   * - ``set_preconditioner(psetup, psolve)``
     - Receive the preconditioner callbacks as Python callables ``psetup()``
       and ``psolve(r, z, tol, lr)``. Either may be ``None``. Implements
       :c:func:`SUNLinSolSetPreconditioner`.
   * - ``set_scaling_vectors(s1, s2)``
     - Receive the left and right scaling ``N_Vector`` objects; either may be
       ``None``. Implements :c:func:`SUNLinSolSetScalingVectors`.
   * - ``set_zero_guess(onoff)``
     - Note whether ``x`` will arrive as the zero vector, allowing the initial
       matrix-vector product to be skipped. Implements
       :c:func:`SUNLinSolSetZeroGuess`.
   * - ``num_iters()``
     - Return the number of iterations used in the last solve. Implements
       :c:func:`SUNLinSolNumIters`.
   * - ``res_norm()``
     - Return the final residual norm from the last solve. Implements
       :c:func:`SUNLinSolResNorm`.
   * - ``resid()``
     - Return the residual ``N_Vector`` from the last solve. Implements
       :c:func:`SUNLinSolResid`.


.. _Python.Usage.CustomObjects.SUNNonlinearSolver:

CustomSUNNonlinearSolver
^^^^^^^^^^^^^^^^^^^^^^^^

``CustomSUNNonlinearSolver(sunctx, solver_type)`` implements the
:c:type:`SUNNonlinearSolver` class. ``solver_type`` is
``SUNNONLINEARSOLVER_ROOTFIND`` for a solver of :math:`F(y) = 0`,
``SUNNONLINEARSOLVER_FIXEDPOINT`` for a solver of :math:`y = G(y)`, or
``SUNNONLINEARSOLVER_HYBRID``. The choice determines which residual convention
the calling package uses, so it must match the iteration implemented in
``solve``.

Only ``solve`` is required.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Method
     - Description
   * - ``solve(y0, y, w, tol, call_lsetup)``
     - **Required.** Solve the nonlinear system, given the predicted value
       ``y0``, the solution ``N_Vector`` ``y`` to fill in, the error weight
       vector ``w``, the convergence tolerance ``tol``, and a flag requesting a
       refresh of the linear solver setup. Implements
       :c:func:`SUNNonlinSolSolve`.
   * - ``initialize()``
     - Perform any setup needed once every callback has been installed.
       Implements :c:func:`SUNNonlinSolInitialize`.
   * - ``setup(y)``
     - Perform per-step setup. Implements :c:func:`SUNNonlinSolSetup`.
   * - ``set_max_iters(maxiters)``
     - Set the maximum number of iterations per solve. Implements
       :c:func:`SUNNonlinSolSetMaxIters`.
   * - ``get_num_iters()``
     - Return ``(status, count)`` for the total iteration count. Implements
       :c:func:`SUNNonlinSolGetNumIters`.
   * - ``get_cur_iter()``
     - Return ``(status, index)`` for the iteration currently in progress -- not
       a total; CVODES and ARKODE read this while deciding whether to refresh
       the linear solver. Implements :c:func:`SUNNonlinSolGetCurIter`.
   * - ``get_num_conv_fails()``
     - Return ``(status, count)`` for the number of convergence failures.
       Implements :c:func:`SUNNonlinSolGetNumConvFails`.
   * - ``set_options(id, file_name, args)``
     - Handle runtime configuration options. Implements
       :c:func:`SUNNonlinSolSetOptions`.

A nonlinear solver is unusual among these interfaces in that it does not
receive the functions it needs as constructor arguments: the calling package
*installs* them by calling setter operations, once, before the first solve. Each
setter receives an ordinary Python callable, so SUNDIALS' function and data
pointers never appear in the subclass's interface. Store the callable and use it
from ``solve``.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Setter
     - Callable it receives
   * - ``set_sys_fn(sys_fn)``
     - ``sys_fn(y, F)`` evaluates the nonlinear residual ``F`` at ``y``.
   * - ``set_sys_fns(root_fn, fixed_point_fn)``
     - Both residual forms at once, for a hybrid solver.
   * - ``set_lsetup_fn(lsetup_fn)``
     - ``lsetup_fn(jbad)`` returns ``(status, jcur)``, forming or refreshing
       the Newton matrix. ``jbad`` reports that the caller believes the current
       matrix is stale; ``jcur`` reports whether the result is freshly computed.
   * - ``set_lsolve_fn(lsolve_fn)``
     - ``lsolve_fn(b)`` solves the Newton linear system in place, overwriting
       ``b`` with the solution.
   * - ``set_conv_test_fn(conv_test_fn)``
     - ``conv_test_fn(y, del, tol, ewt)`` returns ``SUN_SUCCESS`` when
       converged, ``SUN_NLS_CONTINUE`` to keep iterating, or a failure code.
   * - ``set_norm_fn(norm_fn)``
     - ``norm_fn(del, ewt)`` returns ``(status, norm)``.
   * - ``set_get_update_norm_fn(fn)``
     - ``fn()`` returns ``(status, norm)`` for the last update.
   * - ``set_get_conv_rate_fn(fn)``
     - ``fn()`` returns ``(status, rate)`` for the observed convergence rate.

.. warning::

   The callables received from ``set_sys_fn``, ``set_sys_fns``,
   ``set_lsetup_fn``, and ``set_lsolve_fn`` are valid **only while SUNDIALS is
   inside your** ``setup()`` **or** ``solve()``. These four callbacks take an
   opaque memory pointer that the calling package supplies on entry, so outside
   such a call there is no pointer to pass through; invoking one then raises
   ``RuntimeError`` rather than using a stale pointer. The callables from the
   remaining setters carry their own data and may be called at any time.

   Installing a callback also revokes the callable previously installed in that
   slot, matching the semantics of the C setters. A revoked callable raises
   ``RuntimeError`` instead of dispatching through a function pointer SUNDIALS
   has already replaced.


.. _Python.Usage.CustomObjects.SUNAdaptController:

CustomSUNHController and CustomSUNMRIController
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``CustomSUNHController(sunctx)`` and ``CustomSUNMRIController(sunctx)``
implement the :c:type:`SUNAdaptController` class with controller types
``SUN_ADAPTCONTROLLER_H`` and ``SUN_ADAPTCONTROLLER_MRI_H_TOL`` respectively.
Use the first for single-rate time step adaptivity and the second for the
multirate controllers used by MRIStep.

A controller sees very little of the integrator: only the step size just
attempted, the order of the error estimate, and a scaled error measure ``dsm``
formed from the embedded error estimate, where ``dsm`` near 1 means the step met
the requested tolerance, greater than 1 means it failed the error test and will
be retried, and less than 1 means it was more accurate than needed.

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Method
     - Description
   * - ``estimate_step(h, p, dsm)``
     - **Required for** ``CustomSUNHController``. Return ``(status, hnew)``, the
       step size to attempt next. The calling package applies its own step size
       bounds and change ratio limits to the returned value. Implements
       :c:func:`SUNAdaptController_EstimateStep`.
   * - ``estimate_step_tol(H, tolfac, P, DSM, dsm)``
     - **Required for** ``CustomSUNMRIController``. Return
       ``(status, Hnew, tolfacnew)``, the next slow step size and the next fast
       solver tolerance factor. Implements
       :c:func:`SUNAdaptController_EstimateStepTol`.
   * - ``reset()``
     - Discard accumulated history. Called when the integration is
       reinitialized; a controller that keeps history and does not implement
       this will make poor predictions afterwards. Implements
       :c:func:`SUNAdaptController_Reset`.
   * - ``set_defaults()``
     - Restore the controller parameters to their default values. Implements
       :c:func:`SUNAdaptController_SetDefaults`.
   * - ``set_error_bias(bias)``
     - Store the multiplier to apply to ``dsm`` before using it. The controller
       must apply the bias itself. Implements
       :c:func:`SUNAdaptController_SetErrorBias`.
   * - ``update_h(h, dsm)``
     - Record an *accepted* step and its error measure. History belongs here
       rather than in ``estimate_step``, which is also called for rejected
       attempts. Implements :c:func:`SUNAdaptController_UpdateH`.
   * - ``update_mri_h_tol(H, tolfac, DSM, dsm)``
     - The multirate counterpart of ``update_h``. Implements
       :c:func:`SUNAdaptController_UpdateMRIHTol`.
