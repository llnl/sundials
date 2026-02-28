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

.. _SUNDIALS.Logging:

Status and Error Logging
========================

.. versionadded:: 6.2.0

SUNDIALS includes a built-in logging functionality which can be used to direct
error messages, warning messages, informational output, and debugging output to
specified files. This capability requires enabling both build-time and run-time
options to ensure the best possible performance is achieved.

.. _SUNDIALS.Logging.Enabling:

Enabling Logging
----------------

To enable logging, the CMake option :cmakeop:`SUNDIALS_LOGGING_LEVEL` must be
set to the maximum desired output level when configuring SUNDIALS. See the
:cmakeop:`SUNDIALS_LOGGING_LEVEL` documentation for the numeric values
corresponding to errors, warnings, info output, and debug output where errors <
warnings < info output < debug output < extra debug output. By default only
warning and error messages are logged.

When SUNDIALS is built with logging enabled, then the default logger (stored in
the :c:type:`SUNContext` object) may be configured through environment variables
without any changes to user code. The available environment variables are:

.. code-block::

   SUNLOGGER_ERROR_FILENAME
   SUNLOGGER_WARNING_FILENAME
   SUNLOGGER_INFO_FILENAME
   SUNLOGGER_DEBUG_FILENAME

These environment variables may be set to a filename string. There are two
special filenames: ``stdout`` and ``stderr``. These two filenames will result in
output going to the standard output file and standard error file. For example,
consider the CVODE Roberts example, where we can direct the informational output
to the file ``sun.log`` as follows

.. code-block::

   SUNLOGGER_INFO_FILENAME=sun.log ./examples/cvode/serial/cvRoberts_dns

The different environment variables may all be set to the same file, or to
distinct files, or some combination there of. To disable output for one of the
streams, then do not set the environment variable, or set it to an empty string.
If :cmakeop:`SUNDIALS_LOGGING_LEVEL` was set at build-time to a level lower than
the corresponding environment variable, then setting the environment variable
will do nothing. For example, if the logging level is set to ``2`` (errors and
warnings), setting ``SUNLOGGER_INFO_FILENAME`` will do nothing.

Alternatively, the default logger can be accessed with
:c:func:`SUNContext_GetLogger` and configured using the
:ref:`SUNDIALS.Logging.API` or a user may create, configure, and attach a
non-default logger using the :ref:`SUNDIALS.Logging.API`.

.. warning::

   A non-default logger should be created and attached to the context object
   prior to any other SUNDIALS calls in order to capture all log events.

The following examples demonstrate how to use the logging interface via the C
API:

.. code-block::

   examples/arkode/CXX_serial/ark_analytic_sys.cpp
   examples/cvode/serial/cvAdvDiff_bnd.c
   examples/cvode/parallel/cvAdvDiff_diag_p.c
   examples/kinsol/CXX_parallel/kin_em_p.cpp
   examples/kinsol/CUDA_mpi/kin_em_mpicuda.cpp

.. _SUNDIALS.Logging.Output:

Logging Output
--------------

Error or warning logs are a single line output with an error or warning message
of the form

.. code-block:: text

   [level][rank][scope][label] message describing the error or warning

where the values in brackets have the following meaning:

* ``level`` is the log level of the message and will be ``ERROR``, ``WARNING``,
  ``INFO``, or ``DEBUG``

* ``rank`` is the MPI rank the message was written from (``0`` by default or if
  SUNDIALS was built without MPI enabled)

* ``scope`` is the message scope i.e., the name of the function from which the
  message was written

* ``label`` provides additional context or information about the logging
  output e.g., ``begin-step``, ``end-linear-solve``, etc.

Informational or debugging logs are either a single line output with a
comma-separated list of key-value pairs of the form

.. code-block:: text

   [level][rank][scope][label] key1 = value1, key2 = value2

or multiline output with one value per line for keys corresponding to a vector
or array e.g.,

.. code-block:: text

   [level][rank][scope][label] y(:) =
   y[0]
   y[1]
   ...

.. note::

   When extra debugging output is enabled, the output will include vector values
   (so long as the :c:type:`N_Vector` used supports printing). Depending on the
   problem size, this may result in very large logging files.

.. _SUNDIALS.Logging.Tools:

Logging Tools
-------------

.. versionadded:: 7.2.0

To assist with extracting data from logging output files, the ``tools``
directory contains the ``suntools`` Python module which provides utilities for
parsing log files in the ``logs`` sub-module.

.. autofunction:: logs.log_file_to_list

.. autofunction:: logs.print_log

.. autofunction:: logs.get_history

The ``tools`` directory also contains example scripts demonstrating how to use
the log parsing functions to extract and plot data.

* ``log_example_print.py`` -- parses a log file and prints the log file list.

* ``log_example.py`` -- plots the step size, order, or error estimate history
  from an ARKODE, CVODE(S), or IDA(S) log file.

* ``log_example_mri.py`` -- plots the step size history from an ARKODE MRIStep
  log file.

For more details on using an example script, run the script with the ``--help``
flag.

.. _SUNDIALS.Logging.API:

Logger API
----------

The central piece of the Logger API is the :c:type:`SUNLogger` type:

.. c:type:: SUNLogger

   An opaque pointer containing logging information.

When SUNDIALS is built with logging enabled, a default logging object is stored
in the :c:type:`SUNContext` object and can be accessed with a call to
:c:func:`SUNContext_GetLogger`.

The enumerated type :c:enum:`SUNLogLevel` is used by some of the logging
functions to identify the output level or file.

.. c:enum:: SUNLogLevel

   The SUNDIALS logging level

   .. c:enumerator:: SUN_LOGLEVEL_ALL

      Represents all output levels

   .. c:enumerator:: SUN_LOGLEVEL_NONE

      Represents none of the output levels

   .. c:enumerator:: SUN_LOGLEVEL_ERROR

      Represents error-level logging messages

   .. c:enumerator:: SUN_LOGLEVEL_WARNING

      Represents warning-level logging messages

   .. c:enumerator:: SUN_LOGLEVEL_INFO

      Represents info-level logging messages

   .. c:enumerator:: SUN_LOGLEVEL_DEBUG

      Represents deubg-level logging messages


The :c:type:`SUNLogger` class provides the following methods.


.. c:function:: int SUNLogger_Create(SUNComm comm, int output_rank, SUNLogger* logger)

   Creates a new :c:type:`SUNLogger` object.

   **Arguments:**
      * ``comm`` -- the MPI communicator to use, if MPI is enabled, otherwise can be ``SUN_COMM_NULL``.
      * ``output_rank`` -- the MPI rank used for output (can be ``-1`` to print
        to all ranks).
      * ``logger`` -- [in,out] On input this is a pointer to a
        :c:type:`SUNLogger`, on output it will point to a new
        :c:type:`SUNLogger` instance.

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred.


.. c:function:: int SUNLogger_CreateFromEnv(SUNComm comm, SUNLogger* logger)

   Creates a new :c:type:`SUNLogger` object and opens the output streams/files
   from the environment variables:

   .. code-block::

      SUNLOGGER_ERROR_FILENAME
      SUNLOGGER_WARNING_FILENAME
      SUNLOGGER_INFO_FILENAME
      SUNLOGGER_DEBUG_FILENAME

   **Arguments:**
      * ``comm`` -- the MPI communicator to use, if MPI is enabled, otherwise can be   ``SUN_COMM_NULL``.
      * ``logger`` -- [in,out] On input this is a pointer to a
        :c:type:`SUNLogger`, on output it will point to a new
        :c:type:`SUNLogger` instance.

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred.


.. c:function:: int SUNLogger_SetErrorFilename(SUNLogger logger, const char* error_filename)

   Sets the filename for error output.

   **Arguments:**
      * ``logger`` -- a :c:type:`SUNLogger` object.
      * ``error_filename`` -- the name of the file to use for error output.

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred.


.. c:function:: int SUNLogger_SetWarningFilename(SUNLogger logger, const char* warning_filename)

   Sets the filename for warning output.

   **Arguments:**
      * ``logger`` -- a :c:type:`SUNLogger` object.
      * ``warning_filename`` -- the name of the file to use for warning output.

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred.


.. c:function:: int SUNLogger_SetInfoFilename(SUNLogger logger, const char* info_filename)

   Sets the filename for info output.

   **Arguments:**
      * ``logger`` -- a :c:type:`SUNLogger` object.
      * ``info_filename`` -- the name of the file to use for info output.

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred.


.. c:function:: int SUNLogger_SetDebugFilename(SUNLogger logger, const char* debug_filename)

   Sets the filename for debug output.

   **Arguments:**
      * ``logger`` -- a :c:type:`SUNLogger` object.
      * ``debug_filename`` -- the name of the file to use for debug output.

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred.


.. c:function:: int SUNLogger_QueueMsg(SUNLogger logger, SUNLogLevel lvl, const char* scope, const char* label, const char* msg_txt, ...)

   Queues a message to the output log level.

   **Arguments:**
      * ``logger`` -- a :c:type:`SUNLogger` object.
      * ``lvl`` -- the message log level (i.e. error, warning, info, debug).
      * ``scope`` -- the message scope (e.g. the function name).
      * ``label`` -- the message label.
      * ``msg_txt`` -- the message text itself.
      * ``...`` -- the format string arguments

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred.


.. c:function:: int SUNLogger_Flush(SUNLogger logger, SUNLogLevel lvl)

   Flush the message queue(s).

   **Arguments:**
      * ``logger`` -- a :c:type:`SUNLogger` object.
      * ``lvl`` -- the message log level (i.e. error, warning, info, debug or
        all).

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred.


.. c:function:: int SUNLogger_GetOutputRank(SUNLogger logger, int* output_rank)

   Get the output MPI rank for the logger.

   **Arguments:**
      * ``logger`` -- a :c:type:`SUNLogger` object.
      * ``output_rank`` -- [in,out] On input this is a pointer to an int, on
        output it points to the int holding the output rank.

   **Returns:**
      * Returns zero if successful, or non-zero if an error occurred.


.. c:function:: int SUNLogger_Destroy(SUNLogger* logger)

   Free the memory for the :c:type:`SUNLogger` object.

   **Arguments:**
      * ``logger`` -- a pointer to the :c:type:`SUNLogger` object.

   **Returns:**
      * Returns zero if successful, or non-zero if an error occur.
