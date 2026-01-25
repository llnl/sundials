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
using SUNDIALS from C or C++. The few notable differences are discussed below.

View Classes and Memory Management
----------------------------------

sundials4py provides natural usage of SUNDIALS objects with object lifetimes managed by the Python garbage collection as with any other Python object.
There is only one caveat, the SUNDIALS integrator/solver ``void*`` objects are wrapped in "View" classes (behind the scenes) for compatibility with nanobind.
These view objects cannot be implicitly converted to the underlying ``void*``. As such, when calling a function which operates on these ``void*`` objects, one must
extract the ``void*`` "capsule" from the view object by calling the view's ``get`` method.

.. code-block:: python
   
   from sundials4py.core import *
   from sundials4py.cvode import *

   ode_problem = MyODEProblemClass()

   sunctx = SUNContext_Create(SUN_COMM_NULL)
   
   cvode = CVodeCreate(CV_BDF, sunctx)

   # notice we need to call cvode.get()
   status = CVodeInit(cvode.get(), ode_problem.f, T0, y)


Return-by-Pointer Parameters
----------------------------

Functions that return values via pointer arguments in the C API are mapped to
Python functions that return a tuple where the **first element** is the
function's return value (typically an error code) and **subsequent elements**
are the values that would be returned via pointer arguments in C, in the same
order as the C function signature.

**Example 1: Single Return-by-Pointer Value**

C:
   .. code-block:: C

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

      long int nsteps, nfevals, nlinsetups, netfails;
      retval = CVodeGetIntegratorStats(cvode_mem, &nsteps, &nfevals, &nlinsetups, &netfails);
      printf("Steps: %ld, Function evals: %ld, Linear setups: %ld, Error test fails: %ld\n", numsteps, nfevals, nlinsetups, netfails);

Python:
   .. code-block:: python

      retval, nsteps, nfevals, nlinsetups, netfails = CVodeGetIntegratorStats(cvode_mem.get())
      print(f"Steps: {nsteps}, Function evals: {nfevals}, Linear setups: {nlinsetups}, Error test fails: {netfails}")


Arrays
------

``N_Vector`` objects in sundials4py are compatible with numpy's ``ndarray``. Each ``N_Vector`` can work on a numpy arrays without copies, and you can access
and modify the underlying data directly using ``N_VGetArrayPointer``, which returns a numpy ``ndarray`` view of the data.

SUNDIALS matrix types (dense, banded, sparse) are also exposed as Python objects that provide access to their underlying data as numpy arrays (e.g., via ``SUNDenseMatrix_Data``).

Arrays of scalars (e.g., scaling factors passed to ``N_VLinearCombination``) are also represented as numpy arrays.

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

This allows you to use numpy operations for vector and matrix data, and to pass numpy arrays to and from SUNDIALS routines efficiently and without unnecessary copies.


User-Supplied Callback Functions
--------------------------------

SUNDIALS packages and several modules/classes require user-supplied callback functions to define problem-specific behavior, 
such as the right-hand side of an ODE or a nonlinear system function. In sundials4py, you can provide these as standard Python functions or lambdas.

The callback signatures follow the C API. As such, ``N_Vector`` arguments are passed as ``N_Vector`` objects and the underlying ndarray must be extracted in the user code. The only caveat is that return-by-pointer parameters are removed from the signature, and instead become return values (mirroring how return-by-pointer parameters for other functions are handled)

Most callback signatures include a ``void* user_data`` argument. In Python, this argument must be present in the signature, but it should be ignored.

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


.. warning::

   The ``user_data`` argument should always be ``None`` or ``_`` on the Python side. If it is listed otherwise, then it should be ignored to avoid causing catastrophic errors.


Error Codes
-----------

The named ``SUN_ERR_*`` code constants are not available in Python. However, all negative values of ``SUNErrCode``
are still errors, zero is success, and positive values are warnings. As such, users  Users can call ``SUNGetErrMsg``
from Python with the returned ``SUNErrCode`` to get further information about an error.

Examples
--------

Examples can be found in ``examples/python``.
