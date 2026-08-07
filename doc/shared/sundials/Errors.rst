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

.. _SUNDIALS.Errors:

Error Checking
==============

.. versionadded:: 7.0.0

Until version 7.0.0, error reporting and handling was inconsistent throughout SUNDIALS. Starting with version 7.0.0 all of SUNDIALS (the core, implementations of core modules, and
packages) reports error messages through the :c:type:`SUNLogger` API. Furthermore, functions in the
SUNDIALS core API (i.e., ``SUN`` or ``N_V`` functions only) either return a :c:type:`SUNErrCode`, or
(if they don't return a :c:type:`SUNErrCode`) they internally record an error code (if an error
occurs) within the :c:type:`SUNContext` for the execution stream. This "last error" is accessible
via the :c:func:`SUNContext_GetLastError` or :c:func:`SUNContext_PeekLastError` functions.

.. c:type:: int SUNErrCode

Thus, in user code, SUNDIALS core API functions can be checked for errors in one of two ways:

.. code-block:: C

  SUNContext sunctx;
  SUNErrCode sunerr;
  N_Vector v;
  int length;
  sunrealtype dotprod;

  // Every code that uses SUNDIALS must create a SUNContext.
  sunerr = SUNContext_Create(comm, &sunctx)
  if (sunerr) { /* an error occurred, do something */ }

  // Create a SUNDIALS serial vector.
  // Some functions do not return an error code.
  // We have to check for errors in these functions using SUNContext_GetLastError.
  length = 2;
  v = N_VNew_Serial(length, sunctx);
  sunerr = SUNContext_GetLastError(sunctx);
  if (sunerr) {  /* an error occurred, do something */ }

  // If the function returns a SUNErrCode, we can check it directly
  sunerr = N_VLinearCombination(...);
  if (sunerr) {  /* an error occurred, do something */ }

  // Another function that does not return a SUNErrCode.
  dotprod = N_VDotProd(...);
  SUNContext_GetLastError(sunctx);
  if (sunerr) {
   /* an error occurred, do something */
  } else {
    print("dotprod = %.2f\n", dotprod);
  }


The function :c:func:`SUNGetErrMsg` can be used to get a message describing the error code.

.. c:function:: const char* SUNGetErrMsg(SUNErrCode code)

  Returns a message describing the error code.

  :param code: the error code

  :return: a message describing the error code.


.. note::

  It is recommended in most cases that users check for an error after calling SUNDIALS functions.
  However, users concerned with getting the most performance might choose to exclude or limit these checks.


.. warning::

  If a function returns a :c:type:`SUNErrCode` then the return value is the only place the error is available
  i.e., these functions do not store their error code as the "last error" so it is invalid to use
  :c:func:`SUNContext_GetLastError` to check these functions for errors.


.. _SUNDIALS.Errors.Handlers:

Error Handler Functions
-----------------------

When an error occurs in SUNDIALS, it calls error handler functions that have
been pushed onto the error handler stack in last-in first-out order.
Specific error handlers can be enabled by pushing them onto the error handler stack
with the function :c:func:`SUNContext_PushErrHandler`. They may disabled by calling :c:func:`SUNContext_PopErrHandler` or :c:func:`SUNContext_ClearErrHandlers`.
A SUNDIALS error handler function has the type

.. c:type:: void (*SUNErrHandlerFn)(int line, const char* func, const char* file, \
                                           const char* msg, SUNErrCode err_code, \
                                           void* err_user_data, SUNContext sunctx)

SUNDIALS provides a few different error handlers that can be used, or a custom one defined by the
user can be provided (useful for linking SUNDIALS errors to your application's error handling).
The default error handler is :c:func:`SUNLogErrHandlerFn` which logs an error to a specified
file or ``stderr`` if no file is specified.

The error handlers provided in SUNDIALS are:

.. c:function:: void SUNLogErrHandlerFn(int line, const char* func, const char* file, \
                                        const char* msg, SUNErrCode err_code, \
                                        void* err_user_data, SUNContext sunctx)

  Logs the error that occurred using the :c:type:`SUNLogger` from ``sunctx``.
  This is the default error handler.

  :param line: the line number at which the error occurred
  :param func: the function in which the error occurred
  :param file: the file in which the error occurred
  :param msg: the message to log, if this is ``NULL`` then the default error message for the error code will be used
  :param err_code: the error code for the error that occurred
  :param err_user_data: the user pointer provided to :c:func:`SUNContext_PushErrHandler`
  :param sunctx: pointer to a valid :c:type:`SUNContext` object

  :return: ``void``

.. c:function:: void SUNAbortErrHandlerFn(int line, const char* func, const char* file, \
                                          const char* msg, SUNErrCode err_code, \
                                          void* err_user_data, SUNContext sunctx)

  Logs the error and aborts the program if an error occurred.

  :param line: the line number at which the error occurred
  :param func: the function in which the error occurred
  :param file: the file in which the error occurred
  :param msg: this parameter is ignored
  :param err_code: the error code for the error that occurred
  :param err_user_data: the user pointer provided to :c:func:`SUNContext_PushErrHandler`
  :param sunctx: pointer to a valid :c:type:`SUNContext` object

  :return: ``void``


.. c:function:: void SUNMPIAbortErrHandlerFn(int line, const char* func, const char* file, \
                                             const char* msg, SUNErrCode err_code, \
                                             void* err_user_data, SUNContext sunctx)

  Logs the error and calls ``MPI_Abort`` if an error occurred.

  :param line: the line number at which the error occurred
  :param func: the function in which the error occurred
  :param file: the file in which the error occurred
  :param msg: this parameter is ignored
  :param err_code: the error code for the error that occurred
  :param err_user_data: the user pointer provided to :c:func:`SUNContext_PushErrHandler`
  :param sunctx: pointer to a valid :c:type:`SUNContext` object

  :return: ``void``
