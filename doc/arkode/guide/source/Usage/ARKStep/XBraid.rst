.. ----------------------------------------------------------------
   Programmer(s): David J. Gardner @ LLNL
   ----------------------------------------------------------------
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

.. _ARKODE.Usage.ARKStep.XBraid:

Multigrid Reduction in Time with XBraid
=======================================

The prior sections discuss using ARKODE in a traditional sequential
time integration setting i.e., the solution is advanced from one time to the
next where all parallelism resides within the evaluation of a step e.g., the
computation of the right-hand side, (non)linear solves, vector operations etc.
For example, when discretizing a partial differential equation using a method
of lines approach the spatially-discretized equations comprise a large set
of ordinary differential equations that can be evolved with ARKODE. In this
case the parallelization lies in decomposing the spatial domain unknowns across
distributed computational nodes. Considering the strong scaling case at a given
spatial resolution, as the problem is spread across greater numbers of
computational nodes scalability in the spatial dimension is exhausted and
sequential time integration becomes a bottleneck. This bottleneck is largely
driven by the hardware shift from faster clock speeds to greater concurrency to
achieve performance gains. In this case, at the spatial scaling limit and with
stagnant clock speeds, more time steps will lead to an increased runtime.

An alternative approach to sequential time integration is to solve for all time
values simultaneously. One such approach is multigrid reduction in time :cite:p:`FFKMS:14`
(MGRIT) which uses a highly parallel iterative method to expose parallelism in
the time domain in addition to the spatial parallelization. Starting with an
initial temporal grid the multilevel algorithm constructs successively coarser
time grids and uses each coarse grid solution to improve the solution at the
next finer scale. In the two level case the MGRIT algorithm is as follows:

#. Relax the solution on the fine grid (parallel-in-time)

#. Restrict the solution to the fine grid (time re-discretization).

#. Solve the residual equation on the coarse grid (serial-in-time).

#. Correct the fine grid solution (parallel-in-time).

Applying this algorithm recursively for the solve step above leads to the
multilevel algorithm.

The XBraid library :cite:p:`xbraid` implements the MGRIT algorithm in a
non-intrusive manner, enabling the reuse of existing software for sequential
time integration. The following sections describe the ARKODE + XBraid interface
and the steps necessary to modify an existing code that already uses ARKODE's
ARKStep time-stepping module to also use XBraid.



.. _ARKODE.Usage.ARKStep.SUNBraidInterface:

SUNBraid Interface
------------------

Interfacing ARKODE with XBraid requires defining two data structures. The
first is the XBraid application data structure that contains the data necessary
for carrying out a time step and is passed to every interface function (much
like the user data pointer in SUNDIALS packages). For this structure the
SUNBraid interface defines the generic SUNBraidApp structure described below
that serves as the basis for creating integrator-specific or user-defined
interfaces to XBraid. The second structure holds the problem state data at a
certain time value. This structure is defined by the SUNBraidVector structure
and simply contains an N_Vector. In addition to the two data structures several
functions defined by the XBraid API are required. These functions include vector
operations (e.g., computing vector sums or norms) as well as functions to
initialize the problem state, access the current solution, and take a time step.

The ARKBraid interface, built on the SUNBraidApp and SUNBraidVector structures,
provides all the functionality needed combine ARKODE and XBraid for
parallel-in-time integration. As such, only a minimal number of changes are
necessary to update an existing code that uses ARKODE to also use XBraid.



.. _ARKODE.Usage.ARKStep.SUNBraidApp:

SUNBraidApp
^^^^^^^^^^^

As mentioned above the SUNBraid interface defines the SUNBraidApp structure to
hold the data necessary to compute a time step. This structure, like other
SUNDIALS generic objects, is defined as a structure consisting of an
implementation specific *content* field and an operations structure comprised
of a set of function pointers for implmentation-defined operations on the
object. Specifically the SUNBraidApp type is defined as

.. code-block:: C

   /* Define XBraid App structure */
   struct _braid_App_struct
   {
     void        *content;
     SUNBraidOps ops;
   };

   /* Pointer to the interface object (same as braid_App) */
   typedef struct _braid_App_struct *SUNBraidApp;

Here, the SUNBraidOps structure is defined as

.. code-block:: C

   /* Structure containing function pointers to operations */
   struct _SUNBraidOps
   {
     int (*getvectmpl)(braid_App app, N_Vector *tmpl);
   };

   /* Pointer to operations structure */
   typedef struct _SUNBraidOps *SUNBraidOps;

The generic SUNBraidApp defines and implements the generic operations acting on
a SUNBraidApp object. These generic functions are nothing but wrappers to access
the specific implementation through the object's operations structure. To
illustrate this point we show below the implementation of the
:c:func:`SUNBraidApp_GetVecTmpl()` function:

.. code-block:: C

   /* Get a template vector from the integrator */
   int SUNBraidApp_GetVecTmpl(braid_App app, N_Vector *y)
   {
     if (app->ops->getvectmpl == NULL) return SUNBRAID_OPNULL;
     return app->ops->getvectmpl(app, y);
   }

The SUNBraidApp operations are define below in
:numref:`ARKODE.Usage.ARKStep.SUNBraidOps`.



.. _ARKODE.Usage.ARKStep.SUNBraidOps:

SUNBraidOps
^^^^^^^^^^^

In this section we define the SUNBraidApp operations and, for each operation, we
give the function signature, a description of the expected behavior, and an
example usage of the function.

.. c:function:: int SUNBraidApp_GetVecTmpl(braid_App app, N_Vector *y)

   This function returns a vector to use as a template for creating new vectors
   with :c:func:`N_VClone()`.

   :param app: input, a SUNBraidApp instance (XBraid app structure).
   :param y: output, the template vector.

   :return: If this function is not implemented by the SUNBraidApp
            implementation (i.e., the function pointer is ``NULL``) then this function
            will return *SUNBRAID_OPNULL*. Otherwise the return value depends on the
            particular SUNBraidApp implementation. Users are encouraged to utilize the
            return codes  defined in ``sundials/sundials_xbraid.h`` and listed in
            :numref:`ARKODE.Usage.ARKStep.SUNBraidReturnCodes.Table`.

   .. code-block:: C

      /* Get template vector */
      flag = SUNBraidApp_GetVecTmpl(app, y_ptr);
      if (flag != SUNBRAID_SUCCESS) return flag;



.. _ARKODE.Usage.ARKStep.SUNBraidApp_Utilities:

SUNBraidApp Utility Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition to the generic SUNBraidApp operations the following utility
functions are provided to assist in creating and destroying a SUNBraidApp
instance.

.. c:function:: int SUNBraidApp_NewEmpty(braid_App *app)

   This function creates a new SUNBraidApp instance with the content and
   operations initialized to ``NULL``.

   :param app: output, an empty SUNBraidApp instance (XBraid app structure).

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ALLOCFAIL: if a memory allocation failed.

   .. code-block:: C

      /* Create empty XBraid interface object */
      flag = SUNBraidApp_NewEmpty(app_ptr);
      if (flag != SUNBRAID_SUCCESS) return flag;



.. c:function:: int SUNBraidApp_FreeEmpty(braid_App *app)

   This function destroys an empty SUNBraidApp instance.

   :param app: input, an empty SUNBraidApp instance (XBraid app structure).

   :retval SUNBRAID_SUCCESS: if successful.

   .. code-block:: C

      /* Free empty XBraid interface object */
      flag = SUNBraidApp_FreeEmpty(app_ptr);


   .. warning::

      This function does not free the SUNBraidApp object's content structure. An
      implementation should free its content before calling
      :c:func:`SUNBraidApp_FreeEmpty()` to deallocate the base SUNBraidApp
      structure.



.. _ARKODE.Usage.ARKStep.SUNBraidVector:

SUNBraidVector
^^^^^^^^^^^^^^

As mentioned above the SUNBraid interface defines the SUNBraidVector structure
to store a snapshot of solution data at a single point in time and this
structure simply contains an N_Vector. Specifically, the structure is defined
as follows:

.. c:type:: struct _braid_Vector_struct *SUNBraidVector;

   Pointer to vector wrapper (same as braid_Vector)

.. c:struct:: _braid_Vector_struct

   .. c:member:: N_Vector y

      SUNDIALS N_Vector wrapped by the ``braid_Vector``

To assist in creating creating and destroying this structure the following
utility functions are provided.

.. c:function:: int SUNBraidVector_New(N_Vector y, SUNBraidVector *u)

   This function creates a new SUNBraidVector wrapping the N_Vector y.

   :param y: input, the N_Vector to wrap.
   :param u: output, the SUNBraidVector wrapping *y*.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if *y* is ``NULL``.
   :retval SUNBRAID_ALLOCFAIL: if a memory allocation fails.

   .. code-block:: C

      /* Create new vector wrapper */
      flag = SUNBraidVector_New(y, u_ptr);
      if (flag != SUNBRAID_SUCCESS) return flag;

   .. warning::

      The SUNBraidVector takes ownership of the wrapped N_Vector and as such the
      wrapped N_Vector is destroyed when the SUNBraidVector is freed with
      :c:func:`SUNBraidVector_Free()`.



.. c:function:: int SUNBraidVector_GetNVector(SUNBraidVector u, N_Vector *y)

   This function retrieves the wrapped N_Vector from the SUNBraidVector.

   :param u: input, the SUNBraidVector wrapping *y*.
   :param y: output, the wrapped N_Vector.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if *u* is ``NULL``.
   :retval SUNBRAID_MEMFAIL: if *y* is ``NULL``.

   .. code-block:: C

      /* Create new vector wrapper */
      flag = SUNBraidVector_GetNVector(u, y_ptr);
      if (flag != SUNBRAID_SUCCESS) return flag;



Finally, the SUNBraid interface defines the following vector operations acting
on SUNBraidVectors, that consist of thin wrappers to compatible SUNDIALS
N_Vector operations.

.. c:function:: int SUNBraidVector_Clone(braid_App app, braid_Vector u, braid_Vector *v_ptr)

   This function creates a clone of the input SUNBraidVector and copies the
   values of the input vector *u* into the output vector *v_ptr* using
   :c:func:`N_VClone()` and :c:func:`N_VScale()`.

   :param app: input, a SUNBraidApp instance (XBraid app structure).
   :param u: input, the SUNBraidVector to clone.
   :param v_ptr: output, the new SUNBraidVector.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if *u* is ``NULL``.
   :retval SUNBRAID_MEMFAIL: if the N_Vector *y* wrapped by *u* is ``NULL``.
   :retval SUNBRAID_ALLOCFAIL: if a memory allocation fails.



.. c:function:: int SUNBraidVector_Free(braid_App app, braid_Vector u)

   This function destroys the SUNBraidVector and the wrapped N_Vector
   using :c:func:`N_VDestroy()`.

   :param app: input, a SUNBraidApp instance (XBraid app structure).
   :param u: input, the SUNBraidVector to destroy.

   :retval SUNBRAID_SUCCESS: if successful.



.. c:function:: int SUNBraidVector_Sum(braid_App app, braid_Real alpha, braid_Vector x, braid_Real beta, braid_Vector y)

   This function computes the vector sum
   :math:`\alpha x + \beta y \rightarrow y` using :c:func:`N_VLinearSum()`.

   :param app: input, a SUNBraidApp instance (XBraid app structure).
   :param alpha: input, the constant :math:`\alpha`.
   :param x: input, the vector :math:`x`.
   :param beta: input, the constant :math:`\beta`.
   :param y: input/output, the vector :math:`y`.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if *x* or *y* is ``NULL``.
   :retval SUNBRAID_MEMFAIL: if either of the wrapped N_Vectors are ``NULL``.



.. c:function:: int SUNBraidVector_SpatialNorm(braid_App app, braid_Vector u, braid_Real *norm_ptr)

   This function computes the 2-norm of the vector *u* using
   :c:func:`N_VDotProd()`.

   :param app: input, a SUNBraidApp instance (XBraid app structure).
   :param u: input, the vector *u*.
   :param norm_ptr: output, the L2 norm of *u*.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if *u* is ``NULL``.
   :retval SUNBRAID_MEMFAIL: if the wrapped N_Vector is ``NULL``.



.. c:function:: int SUNBraidVector_BufSize(braid_App app, braid_Int *size_ptr, braid_BufferStatus bstatus)

   This function returns the buffer size for messages to exchange vector data
   using :c:func:`SUNBraidApp_GetVecTmpl` and :c:func:`N_VBufSize()`.

   :param app: input, a SUNBraidApp instance (XBraid app structure).
   :param size_ptr: output, the buffer size.
   :param bstatus: input, a status object to query for information on the message
                   type.

   :retval SUNBRAID_SUCCESS: if successful
   :retval otherwise: an error flag from :c:func:`SUNBraidApp_GetVecTmpl`
                      or :c:func:`N_VBufSize()`.



.. c:function:: int SUNBraidVector_BufPack(braid_App app, braid_Vector u, void *buffer, braid_BufferStatus bstatus)

   This function packs the message buffer for exchanging vector data using
   :c:func:`N_VBufPack()`.

   :param app: input, a SUNBraidApp instance (XBraid app structure).
   :param u: input, the vector to pack into the exchange buffer.
   :param buffer: output, the packed exchange buffer to pack.
   :param bstatus: input, a status object to query for information on the message
                   type.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if *u* is ``NULL``.
   :retval otherwise: An error flag from :c:func:`N_VBufPack()`.



.. c:function:: int SUNBraidVector_BufUnpack(braid_App app, void *buffer, braid_Vector *u_ptr, braid_BufferStatus bstatus)

   This function unpacks the message buffer and creates a new N_Vector and
   SUNBraidVector with the buffer data using :c:func:`N_VBufUnpack()`,
   :c:func:`SUNBraidApp_GetVecTmpl`, and :c:func:`N_VClone()`.

   :param app: input, a SUNBraidApp instance (XBraid app structure).
   :param buffer: input, the exchange buffer to unpack.
   :param u_ptr: output, a new SUNBraidVector containing the buffer data.
   :param bstatus: input, a status object to query for information on the message
                   type.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if *buffer* is ``NULL``.
   :retval SUNBRAID_ALLOCFAIL: if a memory allocation fails.
   :retval otherwise: an error flag from :c:func:`SUNBraidApp_GetVecTmpl` and
                      :c:func:`N_VBufUnpack()`.



.. _ARKODE.Usage.ARKStep.SUNBraidReturnCodes:

SUNBraid Return Codes
^^^^^^^^^^^^^^^^^^^^^

The SUNBraid interface return values are given in
:numref:`ARKODE.Usage.ARKStep.SUNBraidReturnCodes.Table`.

.. _ARKODE.Usage.ARKStep.SUNBraidReturnCodes.Table:
.. table:: SUNBraid Return Codes

   +--------------------------+------------+-------------------------------------+
   | Return value name        | Value      | Meaning                             |
   +--------------------------+------------+-------------------------------------+
   | ``SUNBRAID_SUCCESS``     | :math:`0`  | The call/operation was successful.  |
   +--------------------------+------------+-------------------------------------+
   | ``SUNBRAID_ALLOCFAIL``   | :math:`-1` | A memory allocation failed.         |
   +--------------------------+------------+-------------------------------------+
   | ``SUNBRAID_MEMFAIL``     | :math:`-2` | A memory access fail.               |
   +--------------------------+------------+-------------------------------------+
   | ``SUNBRAID_OPNULL``      | :math:`-3` | The SUNBraid operation is ``NULL``. |
   +--------------------------+------------+-------------------------------------+
   | ``SUNBRAID_ILLINPUT``    | :math:`-4` | An invalid input was provided.      |
   +--------------------------+------------+-------------------------------------+
   | ``SUNBRAID_BRAIDFAIL``   | :math:`-5` | An XBraid function failed.          |
   +--------------------------+------------+-------------------------------------+
   | ``SUNBRAID_SUNFAIL``     | :math:`-6` | A SUNDIALS function failed.         |
   +--------------------------+------------+-------------------------------------+



.. _ARKODE.Usage.ARKStep.ARKBraid:

ARKBraid Interface
------------------

This section describes the ARKBraid implementation of a SUNBraidApp for using
the ARKODE's ARKStep time-stepping module with XBraid. The following section
:numref:`ARKODE.Usage.ARKStep.ARKBraid_InitDealloc` describes routines for creating,
initializing, and destroying the ARKODE + XBraid interface, routines for
setting optional inputs, and routines for retrieving data from an ARKBraid
instance. As noted above, interfacing with XBraid requires providing functions
to initialize the problem state, access the current solution, and take a time
step. The default ARKBraid functions for each of these actions are defined in :numref:`ARKODE.Usage.ARKStep.ARKBraid_Interface`  and may be overridden by
user-defined if desired. A skeleton of the user's main or calling program for
using the ARKBraid interface is given in
:numref:`ARKODE.Usage.ARKStep.ARKBraid_Skeleton`. Finally, for advanced users that
wish to create their own SUNBraidApp implementation using ARKODE,
:numref:`ARKODE.Usage.ARKStep.ARKBraid_Utility` describes some helpful
functions available to the user.



.. _ARKODE.Usage.ARKStep.ARKBraid_InitDealloc:

ARKBraid Initialization and Deallocation Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section describes the functions that are called by the user to create,
initialize, and destroy an ARKBraid instance. Each user-callable function
returns ``SUNBRAID_SUCCESS`` (i.e., 0) on a successful call and a negative value
if an error occurred. The possible return codes are given in
:numref:`ARKODE.Usage.ARKStep.SUNBraidReturnCodes.Table`.



.. c:function:: int ARKBraid_Create(void *arkode_mem, braid_App *app)

   This function creates a SUNBraidApp object, sets the content pointer to the
   private ARKBraid interface structure, and attaches the necessary SUNBraidOps
   implementations.

   :param arkode_mem: input, a pointer to an ARKODE memory structure.
   :param app: output, an ARKBraid instance (XBraid app structure).

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: *arkode_mem* is ``NULL``.
   :retval SUNBRAID_ALLOCFAIL: if a memory allocation failed.

   .. warning::

      The ARKBraid interface is ARKStep-specific. Although one could eventually
      construct an XBraid interface to other of ARKODE time-stepping modules
      (e.g., ERKStep or MRIStep), those are not currently supported by this
      implementation.



.. c:function:: int ARKBraid_BraidInit(MPI_Comm comm_w, MPI_Comm comm_t, sunrealtype tstart, sunrealtype tstop, sunindextype ntime, braid_App app, braid_Core *core)

   This function wraps the XBraid ``braid_Init()`` function to create the
   XBraid core memory structure and initializes XBraid with the ARKBraid and
   SUNBraidVector interface functions.

   :param comm_w: input,  the global MPI communicator for space and time.
   :param comm_t: input,  the MPI communicator for the time dimension.
   :param tstart: input,  the initial time value.
   :param tstop: input,  the final time value.
   :param ntime: input,  the initial number of grid points in time.
   :param app: input,  an ARKBraid instance.
   :param core: output, the XBraid core memory structure.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if either MPI communicator is ``MPI_COMM_NULL``,
                              if *ntime* < 2, or if *app* or its content is ``NULL``.
   :retval SUNBRAID_BRAIDFAIL: if the ``braid_Init()`` call fails. The XBraid return
                               value can be retrieved with :c:func:`ARKBraid_GetLastBraidFlag()`.

   .. note::

      If desired, the default functions for vector initialization, accessing the
      solution, taking a time step, and computing the spatial norm should be
      overridden before calling this function.
      See :numref:`ARKODE.Usage.ARKStep.ARKBraid_Set` for more details.

   .. warning::

      The user is responsible for deallocating the XBraid core memory structure
      with the XBraid function ``braid_Destroy()``.



.. c:function:: int ARKBraid_Free(braid_App *app)

   This function deallocates an ARKBraid instance.

   :param app: input, a pointer to an ARKBraid instance.

   :retval SUNBRAID_SUCCESS: if successful.



.. _ARKODE.Usage.ARKStep.ARKBraid_Set:

ARKBraid Set Functions
^^^^^^^^^^^^^^^^^^^^^^

This section describes the functions that are called by the user to set optional
inputs to control the behavior of an ARKBraid instance or to provide alternative
XBraid interface functions. Each user-callable function returns
``SUNBRAID_SUCCESS`` (i.e., 0) on a successful call and a negative value if an
error occurred. The possible return codes are given in
:numref:`ARKODE.Usage.ARKStep.SUNBraidReturnCodes.Table`.



.. c:function:: int ARKBraid_SetStepFn(braid_App app, braid_PtFcnStep step)

   This function sets the step function provided to XBraid (default
   :c:func:`ARKBraid_Step()`).

   :param app: input, an ARKBraid instance.
   :param step: input, an XBraid step function. If *step* is ``NULL``, the
                default function will be used.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if *app* is ``NULL``.
   :retval SUNBRAID_MEMFAIL: if the *app* content is ``NULL``.

   .. note::

      This function must be called prior to :c:func:`ARKBraid_BraidInit()`.



.. c:function:: int ARKBraid_SetInitFn(braid_App app, braid_PtFcnInit init)

   This function sets the vector initialization function provided to XBraid
   (default :c:func:`ARKBraid_Init()`).

   :param app: input, an ARKBraid instance.
   :param init: input, an XBraid vector initialization function. If *init* is
                ``NULL``, the default function will be used.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if *app* is ``NULL``.
   :retval SUNBRAID_MEMFAIL: if the *app* content is ``NULL``.

   .. note::

      This function must be called prior to :c:func:`ARKBraid_BraidInit()`.



.. c:function:: int ARKBraid_SetSpatialNormFn(braid_App app, braid_PtFcnSpatialNorm snorm)

   This function sets the spatial norm function provided to XBraid (default
   :c:func:`SUNBraidVector_SpatialNorm()`).

   :param app: input, an ARKBraid instance.
   :param snorm: input, an XBraid spatial norm function. If *snorm* is ``NULL``,
                 the default function will be used.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if *app* is ``NULL``.
   :retval SUNBRAID_MEMFAIL: if the *app* content is ``NULL``.

   .. note::

      This function must be called prior to :c:func:`ARKBraid_BraidInit()`.



.. c:function:: int ARKBraid_SetAccessFn(braid_App app, braid_PtFcnAccess access)

   This function sets the user access function provided to XBraid (default
   :c:func:`ARKBraid_Access()`).

   :param app: input, an ARKBraid instance.
   :param init: input, an XBraid user access function. If *access* is ``NULL``,
                the default function will be used.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if *app* is ``NULL``.
   :retval SUNBRAID_MEMFAIL: if the *app* content is ``NULL``.

   .. note::

      This function must be called prior to :c:func:`ARKBraid_BraidInit()`.



.. _ARKODE.Usage.ARKStep.ARKBraid_Get:

ARKBraid Get Functions
^^^^^^^^^^^^^^^^^^^^^^

This section describes the functions that are called by the user to retrieve
data from an ARKBraid instance. Each user-callable function returns
``SUNBRAID_SUCCESS`` (i.e., 0) on a successful call and a negative value if an
error occurred. The possible return codes are given in
:numref:`ARKODE.Usage.ARKStep.SUNBraidReturnCodes.Table`.



.. c:function:: int ARKBraid_GetVecTmpl(braid_App app, N_Vector *tmpl)

   This function returns a vector from the ARKODE memory to use as a template
   for creating new vectors with :c:func:`N_VClone()` i.e., this is the ARKBraid
   implementation of :c:func:`SUNBraidApp_GetVecTmpl()`.

   :param app: input, an ARKBraid instance.
   :param tmpl: output, a template vector.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if *app* is ``NULL``.
   :retval SUNBRAID_MEMFAIL: if the *app* content or ARKODE memory is ``NULL``.



.. c:function:: int ARKBraid_GetARKodeMem(braid_App app, void **arkode_mem)

   This function returns the ARKODE memory structure pointer attached with
   :c:func:`ARKBraid_Create()`.

   :param app: input, an ARKBraid instance.
   :param arkode_mem: output, a pointer to the ARKODE memory structure.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if *app* is ``NULL``.
   :retval SUNBRAID_MEMFAIL: if the *app* content or ARKODE memory is ``NULL``.



.. c:function:: int ARKBraid_GetARKStepMem(braid_App app, void **arkode_mem)

   This function returns the ARKStep memory structure pointer attached with
   :c:func:`ARKBraid_Create()`.

   :param app: input, an ARKBraid instance.
   :param arkode_mem: output, a pointer to the ARKStep memory structure.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if *app* is ``NULL``.
   :retval SUNBRAID_MEMFAIL: if the *app* content or ARKStep memory is ``NULL``.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKBraid_GetARKodeMem` instead.


.. c:function:: int ARKBraid_GetUserData(braid_App app, void **user_data)

   This function returns the user data pointer attached with
   :c:func:`ARKodeSetUserData()`.

   :param app: input, an ARKBraid instance.
   :param user_data: output, a pointer to the user data structure.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if *app* is ``NULL``.
   :retval SUNBRAID_MEMFAIL: if the *app* content or ARKODE memory is ``NULL``.



.. c:function:: int ARKBraid_GetLastBraidFlag(braid_App app, int *last_flag)

   This function returns the return value from the most recent XBraid function
   call.

   :param app: input, an ARKBraid instance.
   :param last_flag: output, the XBraid return value.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if *app* is ``NULL``.
   :retval SUNBRAID_MEMFAIL: if the *app* content is ``NULL``.



.. c:function:: int ARKBraid_GetLastARKodeFlag(braid_App app, int *last_flag)

   This function returns the return value from the most recent ARKODE function
   call.

   :param app: input, an ARKBraid instance.
   :param last_flag: output, the ARKODE return value.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if *app* is ``NULL``.
   :retval SUNBRAID_MEMFAIL: if the *app* content is ``NULL``.



.. c:function:: int ARKBraid_GetLastARKStepFlag(braid_App app, int *last_flag)

   This function returns the return value from the most recent ARKStep function
   call.

   :param app: input, an ARKBraid instance.
   :param last_flag: output, the ARKStep return value.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if *app* is ``NULL``.
   :retval SUNBRAID_MEMFAIL: if the *app* content is ``NULL``.

   .. deprecated:: 6.1.0

      Use :c:func:`ARKBraid_GetLastARKodeFlag` instead.


.. c:function:: int ARKBraid_GetSolution(braid_App app, sunrealtype *tout, N_Vector yout)

   This function returns final time and state stored with the default access
   function :c:func:`ARKBraid_Access()`.

   :param app: input, an ARKBraid instance.
   :param last_flag: output, the ARKODE return value.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if *app* is ``NULL``.
   :retval SUNBRAID_MEMFAIL: if the *app* content or the stored vector is ``NULL``.

   .. warning::

      If providing a non-default access function the final time and state are
      not stored within the ARKBraid structure and this function will return an
      error. In this case the user should allocate space to store any desired
      output within the user data pointer attached to ARKODE with
      :c:func:`ARKodeSetUserData()`. This user data pointer can be retrieved
      from the ARKBraid structure with :c:func:`ARKBraid_GetUserData()`.




.. _ARKODE.Usage.ARKStep.ARKBraid_Interface:

ARKBraid Interface Functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section describes the default XBraid interface functions provided by
ARKBraid and called by XBraid to perform certain actions. Any or all of these
functions may be overridden by supplying a user-defined function through the set
functions defined in :numref:`ARKODE.Usage.ARKStep.ARKBraid_Set`. Each default
interface function returns ``SUNBRAID_SUCCESS`` (i.e., 0) on a successful call
and a negative value if an error occurred. The possible return codes are given
in :numref:`ARKODE.Usage.ARKStep.SUNBraidReturnCodes.Table`.



.. c:function:: int ARKBraid_Step(braid_App app, braid_Vector ustop, braid_Vector fstop, braid_Vector u, braid_StepStatus status)

   This is the default step function provided to XBraid. The step function is
   called by XBraid to advance the vector *u* from one time to the next using
   the ARStep memory structure provided to :c:func:`ARKBraid_Create()`. A
   user-defined step function may be set with :c:func:`ARKBraid_SetStepFn()`.

   :param app: input, an ARKBraid instance.
   :param ustop: input, *u* vector at the new time *tstop*.
   :param fstop: input, the right-hand side vector at the new time *tstop*.
   :param u: input/output, on input the vector at the start time and on return the
            vector at the new time.
   :param status: input, a status object to query for information about *u* and
                  to steer XBraid e.g., for temporal refinement.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if *app* is ``NULL``.
   :retval SUNBRAID_MEMFAIL: if the *app* content or ARKODE memory is ``NULL``.
   :retval SUNBRAID_BRAIDFAIL: if an XBraid function fails. The return value can be
                               retrieved with :c:func:`ARKBraid_GetLastBraidFlag()`.
   :retval SUNBRAID_SUNFAIL: if a SUNDIALS function fails. The return value can be
                             retrieved with :c:func:`ARKBraid_GetLastARKStepFlag()`.

   .. note::

      If providing a non-default implementation of the step function the utility
      function :c:func:`ARKBraid_TakeStep()` should be used to advance the input
      vector *u* to the new time.



.. c:function:: int ARKBraid_Init(braid_App app, sunrealtype t, braid_Vector *u_ptr)

   This is the default vector initialization function provided to XBraid. The
   initialization function is called by XBraid to create a new vector and set
   the initial guess for the solution at time :math:`t`. When using this default
   function the initial guess at all time values is the initial condition
   provided to :c:func:`ARKStepCreate`. A user-defined init function may be
   set with :c:func:`ARKBraid_SetInitFn()`.

   :param app: input, an ARKBraid instance.
   :param t: input, the initialization time for the output vector.
   :param u_ptr: output, the new and initialized SUNBraidVector.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if *app* is ``NULL``.
   :retval SUNBRAID_MEMFAIL: if the *app* content or ARKODE memory is ``NULL``.
   :retval SUNBRAID_ALLOCFAIL: if a memory allocation failed.

   .. note::

      If providing a non-default implementation of the vector initialization
      function the utility functions :c:func:`SUNBraidApp_GetVecTmpl()` and
      :c:func:`SUNBraidVector_New()` can be helpful when creating the new vector
      returned by this function.



.. c:function:: int ARKBraid_Access(braid_App app, braid_Vector u, braid_AccessStatus astatus)

   This is the default access function provided to XBraid. The access function
   is called by XBraid to retrieve the current solution. When using this default
   function the final solution time and state are stored within the ARKBraid
   structure. This information can be retrieved with
   :c:func:`ARKBraid_GetSolution()`. A user-defined access function may be
   set with :c:func:`ARKBraid_SetAccessFn()`.

   :param app: input, an ARKBraid instance.
   :param u: input, the vector to be accessed.
   :param status: input, a status object to query for information about *u*.

   :retval SUNBRAID_SUCCESS: if successful.
   :retval SUNBRAID_ILLINPUT: if any of the inputs are ``NULL``.
   :retval SUNBRAID_MEMFAIL: if the *app* content, the wrapped N_Vector, or the
                             ARKODE memory is ``NULL``.
   :retval SUNBRAID_ALLOCFAIL: if allocating storage for the final solution fails.
   :retval SUNBRAID_BRAIDFAIL: if an XBraid function fails. The return value can be
                               retrieved with :c:func:`ARKBraid_GetLastBraidFlag()`.



.. _ARKODE.Usage.ARKStep.ARKBraid_Skeleton:

A skeleton of the user's main program with XBraid
-------------------------------------------------

In addition to the header files required for the integration of the ODE problem
(see the section :numref:`ARKODE.Usage.Headers`), to use the ARKBraid
interface, the user's program must include the header file
``arkode/arkode_xbraid.h`` which declares the needed function prototypes.

The following is a skeleton of the user's main program (or calling program) for
the integration of an ODE IVP using ARKODE's ARKStep time-stepping module with
XBraid for parallel-in-time integration. Most steps are unchanged from the
skeleton program presented in :numref:`ARKODE.Usage.Skeleton`. New or updated
steps are **bold**.

#. **Initialize MPI**

   If parallelizing in space and time split the global communicator into
   communicators for space and time with ``braid_SplitCommworld()``.

#. *Set problem dimensions*

#. *Set vector of initial values*

#. *Create ARKStep object*

#. *Specify integration tolerances*

#. *Create matrix object*

#. *Create linear solver object*

#. *Set linear solver optional inputs*

#. *Attach linear solver module*

#. *Create nonlinear solver object*

#. *Attach nonlinear solver module*

#. *Set nonlinear solver optional inputs*

#. *Set optional inputs*

#. **Create ARKBraid interface**

   Call the constructor :c:func:`ARKBraid_Create()` to create the XBraid app
   structure.

#. **Set optional ARKBraid inputs**

   See :numref:`ARKODE.Usage.ARKStep.ARKBraid_Set` for ARKBraid inputs.

#. **Initialize the ARKBraid interface**

   Call the initialization function :c:func:`ARKBraid_BraidInit()` to create the
   XBraid core memory structure and attach the ARKBraid interface app and
   functions.

#. **Set optional XBraid inputs**

   See the XBraid documentation for available XBraid options.

#. **Evolve the problem**

   Call ``braid_Drive()`` to evolve the problem with MGRIT.

#. **Get optional outputs**

   See :numref:`ARKODE.Usage.ARKStep.ARKBraid_Get` for ARKBraid outputs.

#. *Deallocate memory for solution vector*

#. *Free solver memory*

#. *Free linear solver memory*

#. **Free ARKBraid and XBraid memory**

   Call :c:func:`ARKBraid_Free()` and ``braid_Destroy`` to deallocate the
   ARKBraid interface and and XBraid core memory structures, respectively.

#. *Finalize MPI*




.. _ARKODE.Usage.ARKStep.ARKBraid_Utility:

Advanced ARKBraid Utility Functions
-----------------------------------

This section describes utility functions utilized in the ARKODE + XBraid
interfacing. These functions are used internally by the above ARKBraid interface
functions but are exposed to the user to assist in advanced usage of
ARKODE and XBraid that requires defining a custom SUNBraidApp implementation.



.. c:function:: int ARKBraid_TakeStep(void *arkode_mem, sunrealtype tstart, sunrealtype tstop, N_Vector y, int *ark_flag)

   This function advances the vector *y* from *tstart* to *tstop* using a
   single ARKODE time step with step size *h = tstop - start*.

   :param arkode_mem: input, the ARKODE memory structure pointer.
   :param tstart: input, the step start time.
   :param tstop: input, the step stop time.
   :param y: input/output, on input the solution a *tstop* and on return, the
             solution at time *tstop* if the step was successful (*ark_flag*
             :math:`\geq 0`) or the solution at time *tstart* if the step failed
             (*ark_flag* < 0).
   :param ark_flag: output, the step status flag. If *ark_flag* is:

                    :math:`= 0` then the step succeeded and, if applicable, met the
                    requested temporal accuracy.

                    :math:`> 0` then the step succeeded but failed to meet the requested
                    temporal accuracy.

                    :math:`< 0` then the step failed e.g., a solver failure occurred.

   :return: If all ARKODE function calls are successful the return
            value is *ARK_SUCCESS*, otherwise the return value is the error flag
            returned from the function that failed.
