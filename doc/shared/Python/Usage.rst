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
   particular, third-party libraries are not yet supported.

.. _Python.Usage.Example:

Example Usage
-------------

We now consider a simple CVODE example to illustrate using sundials4py and
highlight some of the differences to using SUNDIALS from C/C++. The items
highlighted below similarly apply to using other SUNDIALS packages. For more
information on usage differences, continue to the :ref:`next section
<Python.Usage.Differences>`. Additional examples can be found in the
``examples/python`` directory of the :examples:`SUNDIALS GitHub repository
<python>`.

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

``N_Vector`` objects in sundials4py are compatible with numpy's
``ndarray``. Each ``N_Vector`` can work on a numpy arrays without copies, and
you can access and modify the underlying data directly using
``N_VGetArrayPointer``, which returns a numpy ``ndarray`` view of the data.

SUNDIALS matrix types (dense, banded, sparse) are also exposed as Python objects
that provide access to their underlying data as numpy arrays (e.g., via
``SUNDenseMatrix_Data``).

Arrays of scalars (e.g., scaling factors passed to ``N_VLinearCombination``) are
also represented as numpy arrays.

**Example: Accessing and modifying an N_Vector**

.. code-block:: python

   y_nvec = N_VNew_Serial(10, sunctx)
   y = N_VGetArrayPointer(y_nvec)
   y[:] = np.linspace(0, 1, 10)  # Set values using numpy

**Example: Using a matrix as a numpy array**

.. code-block:: python

   mat = SUNDenseMatrix(3, 3, sunctx)
   arr = SUNDenseMatrix_Data(mat)
   arr[:] = np.eye(3)  # Set to identity matrix

This allows you to use numpy operations for vector and matrix data, and to pass
numpy arrays to and from SUNDIALS routines efficiently and without unnecessary
copies.


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
      y = N_VGetArrayPointer(y_nvector)
      ydot = N_VGetArrayPointer(ydot_nvector)
      ydot[:] = -y
      return 0

   ark = ARKStepCreate(rhs, None, t0, y, sunctx)

**Example: Nonlinear system for KINSOL**

.. code-block:: python

   # The C signature is:
   # int(N_Vector u, N_Vector g, void* user_data)
   def fp_function(u_nvector, g_nvector, _): # note _ in place of user_data
      # Compute g = F(u)
      u = N_VGetArrayPointer(u_nvector)
      g = N_VGetArrayPointer(g_nvector)
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
