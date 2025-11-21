.. ----------------------------------------------------------------
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
   ----------------------------------------------------------------

.. _SUNDIALS.SUNContext:

The SUNContext Type
===================

.. versionadded:: 6.0.0

All of the SUNDIALS objects (vectors, linear and nonlinear solvers, matrices,
etc.) that collectively form a SUNDIALS simulation, hold a reference to a common
simulation context object defined by the :c:type:`SUNContext` class.

.. c:type:: SUNContext

   An opaque pointer used by SUNDIALS objects for error handling, logging,
   profiling, etc.

Users should create a :c:type:`SUNContext` object prior to any other calls to
SUNDIALS library functions by calling:

.. c:function:: SUNErrCode SUNContext_Create(SUNComm comm, SUNContext* sunctx)

   Creates a :c:type:`SUNContext` object associated with the thread of execution.
   The data of the :c:type:`SUNContext` class is private.

   :param comm: the MPI communicator or ``SUN_COMM_NULL`` if not using MPI.
   :param sunctx: [in,out] upon successful exit, a pointer to the newly created :c:type:`SUNContext` object.

   :return: :c:type:`SUNErrCode` indicating success or failure.

The created :c:type:`SUNContext` object should be provided to the constructor
routines for different SUNDIALS classes/modules e.g.,

.. code-block:: C

   SUNContext sunctx;
   void* package_mem;
   N_Vector x;

   SUNContext_Create(SUN_COMM_NULL, &sunctx);

   package_mem = CVodeCreate(..., sunctx);
   package_mem = IDACreate(..., sunctx);
   package_mem = KINCreate(..., sunctx);
   package_mem = ARKStepCreate(..., sunctx);

   x = N_VNew_<SomeVector>(..., sunctx);

After all other SUNDIALS code, the :c:type:`SUNContext` object should be freed
with a call to:

.. c:function:: SUNErrCode SUNContext_Free(SUNContext* sunctx)

   Frees the :c:type:`SUNContext` object.

   :param sunctx: pointer to a valid :c:type:`SUNContext` object, ``NULL`` upon successful return.

   :return: :c:type:`SUNErrCode` indicating success or failure.

   .. warning::

      When MPI is being used, the :c:func:`SUNContext_Free` must be called prior
      to ``MPI_Finalize``.


The :c:type:`SUNContext` API further consists of the following functions:

.. c:function:: SUNErrCode SUNContext_GetLastError(SUNContext sunctx)

   Gets the last error code set by a SUNDIALS function call. The function
   then resets the last error code to `SUN_SUCCESS`.

   :param sunctx: a valid :c:type:`SUNContext` object.

   :return: the last :c:type:`SUNErrCode` recorded.


.. c:function:: SUNErrCode SUNContext_PeekLastError(SUNContext sunctx)

   Gets the last error code set by a SUNDIALS function call. The function
   *does not* reset the last error code to `SUN_SUCCESS`.

   :param sunctx: a valid :c:type:`SUNContext` object.

   :return: the last :c:type:`SUNErrCode` recorded.


.. c:function:: SUNErrCode SUNContext_PushErrHandler(SUNContext sunctx, SUNErrHandlerFn err_fn, void* err_user_data)

   Pushes a new :c:type:`SUNErrHandlerFn` onto the error handler stack so that it is called when an
   error occurs inside of SUNDIALS.

   :param sunctx: a valid :c:type:`SUNContext` object.
   :param err_fn: a callback function of type :c:type:`SUNErrHandlerFn` to be pushed onto the error handler stack.
   :param err_user_data: a pointer that will be passed back to the callback function when it is called.

   :return: :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNContext_PopErrHandler(SUNContext sunctx)

   Pops the last :c:type:`SUNErrHandlerFn` off of the error handler stack.

   :param sunctx: a valid :c:type:`SUNContext` object.

   :return: :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNContext_ClearErrHandlers(SUNContext sunctx)

   Clears the entire error handler stack. After doing this it is important to push an error handler
   onto the stack with :c:type:`SUNContext_PushErrHandler` otherwise errors will be ignored.

   :param sunctx: a valid :c:type:`SUNContext` object.

   :return: :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNContext_GetProfiler(SUNContext sunctx, SUNProfiler* profiler)

   Gets the :c:type:`SUNProfiler` object associated with the
   :c:type:`SUNContext` object.

   :param sunctx: a valid :c:type:`SUNContext` object.
   :param profiler: [in,out] a pointer to the :c:type:`SUNProfiler` object
        associated with this context; will be ``NULL`` if profiling is not
        enabled.

   :return: :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNContext_SetProfiler(SUNContext sunctx, SUNProfiler profiler)

   Sets the :c:type:`SUNProfiler` object associated with the
   :c:type:`SUNContext` object.

   :param sunctx: a valid :c:type:`SUNContext` object.
   :param profiler: a :c:type:`SUNProfiler` object to associate with this
        context; this is ignored if profiling is not enabled.

   :return: :c:type:`SUNErrCode` indicating success or failure.


.. c:function:: SUNErrCode SUNContext_SetLogger(SUNContext sunctx, SUNLogger logger)

   Sets the :c:type:`SUNLogger` object associated with the :c:type:`SUNContext`
   object.

   :param sunctx: a valid :c:type:`SUNContext` object.
   :param logger: a :c:type:`SUNLogger` object to associate with this
        context; this is ignored if logging is not enabled.

   :return: :c:type:`SUNErrCode` indicating success or failure.

   .. versionadded:: 6.2.0


.. c:function:: SUNErrCode SUNContext_GetLogger(SUNContext sunctx, SUNLogger* logger)

   Gets the :c:type:`SUNLogger` object associated with the :c:type:`SUNContext` object.

   :param sunctx: a valid :c:type:`SUNContext` object.
   :param logger: [in,out] a pointer to the :c:type:`SUNLogger` object associated with this context; will be ``NULL`` if logging is not enabled.  

   :return: :c:type:`SUNErrCode` indicating success or failure.

   .. versionadded:: 6.2.0


.. _SUNDIALS.SUNContext.Threads:

Implications for task-based programming and multi-threading
------------------------------------------------------------

Applications that need to have *concurrently initialized* SUNDIALS simulations
need to take care to understand the following:

#. A :c:type:`SUNContext` object must only be associated with *one* SUNDIALS simulation
   (a solver object and its associated vectors etc.) at a time.

   - Concurrently initialized is not the same as concurrently executing. Even if
     two SUNDIALS simulations execute sequentially, if both are initialized
     at the same time with the same :c:type:`SUNContext`, behavior is undefined.

   - It is OK to reuse a :c:type:`SUNContext` object with another SUNDIALS simulation
     after the first simulation has completed and all of the simulation's
     associated objects (vectors, matrices, algebraic solvers, etc.) have been
     destroyed.

#. The creation and destruction of a :c:type:`SUNContext` object is cheap, especially
   in comparison to the cost of creating/destroying a SUNDIALS solver object.

The following (incomplete) code examples demonstrate these points using CVODE as
the example SUNDIALS package.

.. code-block:: c

   SUNContext sunctxs[num_threads];
   int cvode_initialized[num_threads];
   void* cvode_mem[num_threads];

   // Create
   for (int i = 0; i < num_threads; i++) {
      sunctxs[i] = SUNContext_Create(...);
      cvode_mem[i] = CVodeCreate(..., sunctxs[i]);
      cvode_initialized[i] = 0; // not yet initialized
      // set optional cvode inputs...
   }

   // Solve
   #pragma omp parallel for
   for (int i = 0; i < num_problems; i++) {
      int retval = 0;
      int tid = omp_get_thread_num();
      if (!cvode_initialized[tid]) {
         retval = CVodeInit(cvode_mem[tid], ...);
         cvode_initialized[tid] = 1;
      } else {
         retval = CVodeReInit(cvode_mem[tid], ...);
      }
      CVode(cvode_mem[i], ...);
   }

   // Destroy
   for (int i = 0; i < num_threads; i++) {
      // get optional cvode outputs...
      CVodeFree(&cvode_mem[i]);
      SUNContext_Free(&sunctxs[i]);
   }

Since each thread has its own unique CVODE and SUNContext object pair, there
should be no thread-safety issues. Users should be sure that you apply the same
idea to the other SUNDIALS objects needed as well (e.g. an ``N_Vector``).

The variation of the above code example demonstrates another possible approach:

.. code-block:: c

   // Create, Solve, Destroy
   #pragma omp parallel for
   for (int i = 0; i < num_problems; i++) {
      int retval = 0;
      void* cvode_mem;
      SUNContext sunctx;

      sunctx = SUNContext_Create(...);
      cvode_mem = CVodeCreate(..., sunctx);
      retval = CVodeInit(cvode_mem, ...);

      // set optional cvode inputs...

      CVode(cvode_mem, ...);

      // get optional cvode outputs...

      CVodeFree(&cvode_mem);
      SUNContext_Free(&sunctx);
   }

So long as the overhead of creating/destroying the CVODE object is small
compared to the cost of solving the ODE, this approach is a fine alternative to
the first approach since :c:func:`SUNContext_Create` and
:c:func:`SUNContext_Free` are much cheaper than the CVODE create/free routines.


.. _SUNDIALS.SUNContext.CPP:

Convenience class for C++ Users
-------------------------------

For C++ users a RAII safe class, ``sundials::Context``, is provided:

.. code-block:: cpp

   namespace sundials {

   class Context : public sundials::ConvertibleTo<SUNContext>
   {
   public:
   explicit Context(SUNComm comm = SUN_COMM_NULL)
   {
      sunctx_ = std::make_unique<SUNContext>();
      SUNContext_Create(comm, sunctx_.get());
   }

   /* disallow copy, but allow move construction */
   Context(const Context&) = delete;
   Context(Context&&)      = default;

   /* disallow copy, but allow move operators */
   Context& operator=(const Context&) = delete;
   Context& operator=(Context&&) = default;

   SUNContext get() override
   {
      return *sunctx_.get();
   }
   SUNContext get() const override
   {
      return *sunctx_.get();
   }
   operator SUNContext() override
   {
      return *sunctx_.get();
   }
   operator SUNContext() const override
   {
      return *sunctx_.get();
   }

   ~Context()
   {
      if (sunctx_) SUNContext_Free(sunctx_.get());
   }

   private:
   std::unique_ptr<SUNContext> sunctx_;
   };

   } // namespace sundials
