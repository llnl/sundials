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

.. _SourceCode.Style:

Style
=====

In this section we describe the style conventions and guidelines for SUNDIALS
source code.

Formatting
----------

All new code added to SUNDIALS should be formatted with `clang-format
<https://clang.llvm.org/docs/ClangFormat.html>`_ for C/C++, `fprettify
<https://github.com/fortran-lang/fprettify>`_ for Fortran, `cmake-format
<https://cmake-format.readthedocs.io>`_ for CMake, and `black
<https://black.readthedocs.io>`_ for Python. The ``.clang-format`` file in the
root of the project defines our configuration for clang-format. We use the
default fprettify settings, except we use 2-space indentation. The
``.cmake-format.py`` file in the root of the project defines our configuration
for cmake-format. We also use the default black settings.


To apply ``clang-format``, ``fprettify``, ``cmake-format``, and ``black`` you
can run:

.. code-block:: shell

   ./scripts/format.sh <path to directories or files to format>

.. warning::

   The output of ``clang-format`` is sensitive to the ``clang-format`` version. We recommend
   that you use version ``17.0.4``, which can be installed from source or with Spack. Alternatively,
   when you open a pull request on GitHub, an action will run ``clang-format`` on the code. If any
   formatting is required, the action will fail. Commenting with the magic keyword ``/autofix`` will
   kick off a GitHub action which will automatically apply the formatting changes needed.

If clang-format breaks lines in a way that is unreadable, use ``//`` to break the line. For example,
sometimes (mostly in C++ code) you may have code like this:

.. code-block:: cpp

   MyClass::callAFunctionOfSorts::doSomething().doAnotherThing().doSomethingElse();

That you would like to format as (for readability):

.. code-block:: cpp

   MyObject::callAFunctionOfSorts()
         .doSomething()
         .doAnotherThing()
         .doSomethingElse();

Clang-format might produce something like:

.. code-block:: cpp

   MyObject::callAFunctionOfSorts().doSomething().doAnotherThing()
         .doSomethingElse();


unless you add the ``//``

.. code-block:: cpp

   MyObject::callAFunctionOfSorts()
         .doSomething()       //
         .doAnotherThing()    //
         .doSomethingElse();  //

There are other scenarios (e.g., a function call with a lot of parameters) where
doing this type of line break is useful for readability too.

.. It may be necessary to override clang-tidy at times. This can be done with the
.. ``NOLINT`` magic comments e.g.,

.. .. code-block:: cpp

..   template<class GkoSolverType, class GkoMatrixType>
..   int SUNLinSolFree_Ginkgo(SUNLinearSolver S)
..   {
..     auto solver{static_cast<LinearSolver<GkoSolverType, GkoMatrixType>*>(S->content)};
..     delete solver; // NOLINT
..     return SUNLS_SUCCESS;
..   }

..   class BaseObject {
..   protected:
..     // NOLINTNEXTLINE(cppcoreguidelines-non-private-member-variables-in-classes)
..     SUNContext sunctx_{};
..   };

.. See the clang-tidy documentation for more details.


.. _Style.Output:

Output
------

For consistent formatting of :c:type:`sunrealtype`, the following macros are
available.

.. c:macro:: SUN_FORMAT_E

   A format specifier for scientific notation. This should be used when
   displaying arrays, matrices, and tables where fixed width alignment aids with
   readability.

   **Example usage:**

   .. code-block:: C

      for (i = 0; i < N; i++) {
         fprintf(outfile, SUN_FORMAT_E "\n", xd[i]);
      }

.. c:macro:: SUN_FORMAT_G

   A format specifier for scientific or standard notation, whichever is more
   compact. It is more reader-friendly than :c:macro:`SUN_FORMAT_E` and should
   be used in all cases not covered by that macro.

   **Example usage:**

   .. code-block:: C

      SUNLogInfo(sunctx->logger, "label", "x = " SUN_FORMAT_G, x);

.. c:macro:: SUN_FORMAT_SG

   Like :c:macro:`SUN_FORMAT_G` but with a leading plus or minus sign.


To aid in printing statistics in functions like :c:func:`CVodePrintAllStats`,
the following utility functions are available.

.. c:function:: void sunfprintf_real(FILE* fp, SUNOutputFormat fmt, sunbooleantype start, const char* name, sunrealtype value)

   Writes a :c:type:`sunrealtype` value to a file pointer using the specified
   format.

   :param fp: The output file pointer.
   :param fmt: The output format.
   :param start: :c:macro:`SUNTRUE` if the value is the first in a series of
                 statistics, and :c:macro:`SUNFALSE` otherwise.
   :param name: The name of the statistic.
   :param value: The value of the statistic.

.. c:function:: void sunfprintf_long(FILE* fp, SUNOutputFormat fmt, sunbooleantype start, const char* name, long value)

   Writes a long value to a file pointer using the specified format.

   :param fp: The output file pointer.
   :param fmt: The output format.
   :param start: :c:macro:`SUNTRUE` if the value is the first in a series of
                 statistics, and :c:macro:`SUNFALSE` otherwise.
   :param name: The name of the statistic.
   :param value: The value of the statistic.

.. c:function:: void sunfprintf_long_array(FILE* fp, SUNOutputFormat fmt, sunbooleantype start, const char* name, long* value, size_t count)

   Writes an array of long values to a file pointer using the specified format.

   :param fp: The output file pointer.
   :param fmt: The output format.
   :param start: :c:macro:`SUNTRUE` if the value is the first in a series of
                 statistics, and :c:macro:`SUNFALSE` otherwise.
   :param name: The name of the statistic.
   :param value: Pointer to the array.
   :param count: The number of elements in the array.

.. _Style.Logging:

Logging
-------

Use the :ref:`macros below <Style.Logging.Macros>` to add informational and
debugging messages to SUNDIALS code rather than adding ``#ifdef
SUNDIALS_LOGGING_<level>`` / ``#endif`` blocks containing calls to
:c:func:`SUNLogger_QueueMsg`. Error and warning messages are handled through
package-specific ``ProcessError`` functions or the ``SUNAssert`` and
``SUNCheck`` macros.

The logging macros help ensure messages follow the required format (see
:numref:`SUNDIALS.Logging.Output`) used by the ``suntools`` Python module to
parse logging output (see :numref:`SUNDIALS.Logging.Tools`). For informational
and debugging output the log message payload (the part after the brackets) must
be either be a comma-separated list of key-value pairs with the key and value
separated by an equals sign with a space on either side e.g.,

.. code-block:: C

   /* log an informational message */
   SUNLogInfo(sunctx->logger, "begin-step", "t = " SUN_FORMAT_G ", h = " SUN_FORMAT_G, t, h);

   /* log a debugging message */
   SUNLogDebug(sunctx->logger, "error-estimates",
               "eqm1 = " SUN_FORMAT_G ", eq = " SUN_FORMAT_G ", eqp1 = " SUN_FORMAT_G,
               eqm1, eq, eqp1);

or the name of a vector/array followed by ``(:) =`` with each vector/array entry
written to a separate line e.g., a vector may be logged with

.. code-block:: C

   SUNLogExtraDebugVec(sunctx->logger, "new-solution", ynew, "ynew(:) =");

where the message can contain format specifiers e.g., if ``Fe`` is an array of
vectors you may use

.. code-block:: C

   SUNLogExtraDebugVec(sunctx->logger, "stage-explicit-rhs", Fe[i], "Fe_%d(:) =", i);

The log parser organizes the logging output into Python dictionaries and lists
of dictionaries. To denote different logging regions, ``begin-`` and ``end-``
message labels of the form ``<begin|end>-<region-name>[-list]`` are used. The
log parsing tool will automatically open/close the dictionary (default) or list
of dictionaries (as indicated by the optional ``-list`` suffix) for the region
corresponding to ``region-name``. This convention enables adding new logging
regions to the code without needing to update the logging parser. There are also
three special region markers that are handled separately:

* ``begin|end-step-attempt`` -- opens/closes step attempt dictionaries and
  opens/closes lists within the current step to hold logging output for
  different time levels (i.e., with MRI methods) or problem partitions (i.e.,
  with splitting methods) based on the current time level and partition
  counters.

* ``begin|end-fast-steps`` -- increments/decrements the current time level
  counter.

* ``begin|end-partition-list`` -- increments/decrements the current partition
  counter and opens/closes the ``partitions`` list.

Logging messages that do not start with ``begin-`` or ``end-`` correspond to
output for the active region and are added to the open output dictionary.

.. _Style.Logging.Example:

Logging Example
^^^^^^^^^^^^^^^

As an example consider the information logging output from a CVODE time
step using a fixed-point iteration for the nonlinear solve. The logging line

.. code-block:: C

   SUNLogInfo(CV_LOGGER, "begin-step-attempt", "step = %li, tn = " SUN_FORMAT_G ",
              h = " SUN_FORMAT_G ", q = %d", cv_mem->cv_nst + 1, cv_mem->cv_tn,
              cv_mem->cv_h, cv_mem->cv_q);

creates a new step attempt dictionary and adds the data in the payload to the
dictionary:

.. code-block:: text

   {
     step : 1
     tn : 0.0
     h : 0.002
   }

Inside the step attempt, the logging line

.. code-block:: C

   SUNLogInfo(CV_LOGGER, "begin-nonlinear-solve", "tol = " SUN_FORMAT_G,
              cv_mem->cv_tq[4]);

adds the ``nonlinear-solve`` key to the step attempt dictionary, activates
the nonlinear solver output dictionary, and adds the data in the payload to the
dictionary:

.. code-block:: text

   {
     step : 1
     tn : 0.0
     h : 0.002
     nonlinear-solve :
     {
       tol : 0.01
     }
   }

In the nonlinear solver, the logging line

.. code-block:: C

    SUNLogInfo(NLS->sunctx->logger, "begin-iterations-list", "");

adds the ``iterations`` key to the ``nonlinear-solve`` dictionary and activates
a new output dictionary in the iterations list:

.. code-block:: text

   {
     step : 1
     tn : 0.0
     h : 0.002
     nonlinear-solve :
     {
       tol : 0.01
       iterations :
       [
         { }
       ]
     }
   }

Within the nonlinear solver iteration loop, the logging line

.. code-block:: C

   SUNLogInfo(NLS->sunctx->logger, "nonlinear-iterate",
              "cur-iter = %d, update-norm = " SUN_FORMAT_G,
              FP_CONTENT(NLS)->niters, N_VWrmsNorm(delta, w));

adds information to the active iterations dictionary:

.. code-block:: text

   {
     step : 1
     tn : 0.0
     h : 0.002
     nonlinear-solve :
     {
       tol : 0.01
       iterations :
       [
         {
           cur-iter : 1
           update-norm : 0.02
         }
       ]
     }
   }

At the end of each nonlinear iteration, logging output such as

.. code-block:: C

   SUNLogInfoIf(FP_CONTENT(NLS)->curiter < FP_CONTENT(NLS)->maxiters - 1,
                NLS->sunctx->logger, "end-iterations-list",
                "status = continue");

if the iteration should continue or

.. code-block:: C

   SUNLogInfo(NLS->sunctx->logger, "end-iterations-list", "status = success");

if the solve was successful will close the active output dictionary in the
iterations list and reactivate the output dictionary for the enclosing scope
(i.e., the nonlinear solve region). If additional iterations are necessary, the
logging line

.. code-block:: C

    SUNLogInfo(NLS->sunctx->logger, "begin-iterations-list", "");

will append and activate a new output dictionary to the iterations list. This
continues until the nonlinear solve is complete and a logging line such as

.. code-block:: C

   SUNLogInfo(CV_LOGGER, "end-nonlinear-solve", "status = success, iters = %li",
              nni_inc);

updates and closes the nonlinear solver output dictionary which reactivates the
step attempt dictionary:

.. code-block:: text

   {
     step : 1
     tn : 0.0
     h : 0.002
     nonlinear-solve :
     {
       tol : 0.01
       iterations :
       [
         {
           cur-iter : 1
           update-norm : 0.02
           status : continue
         }
         {
           cur-iter : 2
           update-norm : 0.002
           status : success
         }
       ]
       status : success
       iters : 2
     }
   }

Finally, a logging call such as

.. code-block:: C

   SUNLogInfo(CV_LOGGER, "end-step-attempt",
              "status = success, dsm = " SUN_FORMAT_G, dsm);

closes the step attempt dictionary:

.. code-block:: text

   {
     step : 1
     tn : 0.0
     h : 0.002
     nonlinear-solve :
     {
       tol : 0.01
       iterations :
       [
         {
           cur-iter : 1
           update-norm : 0.02
           status : continue
         }
         {
           cur-iter : 2
           update-norm : 0.002
           status : success
         }
       ]
       status : success
       iters : 2
     }
     status : success
     dsm : 1.0e-3
   }

This processes then repeats for the next step attempt.

.. _Style.Logging.Macros:

Logging Macros
^^^^^^^^^^^^^^

.. versionadded:: 7.2.0

To log informational messages use the following macros:

.. c:macro:: SUNLogInfo(logger, label, msg_txt, ...)

   When information logging is enabled this macro expands to a call to
   :c:func:`SUNLogger_QueueMsg` to log an informational message. Otherwise, this
   expands to nothing.

   :param logger: the :c:type:`SUNLogger` to handle the message.
   :param label: the ``const char*`` message label.
   :param msg_txt: the ``const char*`` message text, may contain format
                   specifiers.
   :param ...: the arguments for format specifiers in ``msg_txt``.

.. c:macro:: SUNLogInfoIf(condition, logger, label, msg_txt, ...)

   When information logging is enabled this macro expands to a conditional call
   to :c:func:`SUNLogger_QueueMsg` to log an informational message. Otherwise,
   this expands to nothing.

   :param condition: a boolean expression that determines if the log message
                     should be queued.
   :param logger: the :c:type:`SUNLogger` to handle the message.
   :param label: the ``const char*`` message label.
   :param msg_txt: the ``const char*`` message text, may contain format.
                   specifiers.
   :param ...: the arguments for format specifiers in ``msg_txt``.

To log debugging messages use the following macros:

.. c:macro:: SUNLogDebug(logger, label, msg_txt, ...)

   When debugging logging is enabled this macro expands to a call to
   :c:func:`SUNLogger_QueueMsg` to log a debug message. Otherwise, this expands
   to nothing.

   :param logger: the :c:type:`SUNLogger` to handle the message.
   :param label: the ``const char*`` message label.
   :param msg_txt: the ``const char*`` message text, may contain format.
                   specifiers.
   :param ...: the arguments for format specifiers in ``msg_txt``.

.. c:macro:: SUNLogDebugIf(condition, logger, label, msg_txt, ...)

   When debugging logging is enabled this macro expands to a conditional call to
   :c:func:`SUNLogger_QueueMsg` to log a debug message. Otherwise, this expands
   to nothing.

   :param condition: a boolean expression that determines if the log message
                     should be queued.
   :param logger: the :c:type:`SUNLogger` to handle the message.
   :param label: the ``const char*`` message label.
   :param msg_txt: the ``const char*`` message text, may contain format.
                   specifiers.
   :param ...: the arguments for format specifiers in ``msg_txt``.

To log extra debugging messages use the following macros:

.. c:macro:: SUNLogExtraDebug(logger, label, msg_txt, ...)

   When extra debugging logging is enabled, this macro expands to a call to
   :c:func:`SUNLogger_QueueMsg` to log an extra debug message. Otherwise, this expands
   to nothing.

   :param logger: the :c:type:`SUNLogger` to handle the message.
   :param label: the ``const char*`` message label.
   :param msg_txt: the ``const char*`` message text, may contain format
                   specifiers.
   :param ...: the arguments for format specifiers in ``msg_txt``.

.. c:macro:: SUNLogExtraDebugIf(condition, logger, label, msg_txt, ...)

   When extra debugging logging is enabled, this macro expands to a conditional
   call to :c:func:`SUNLogger_QueueMsg` to log an extra debug message. Otherwise, this
   expands to nothing.

   :param condition: a boolean expression that determines if the log message
                     should be queued.
   :param logger: the :c:type:`SUNLogger` to handle the message.
   :param label: the ``const char*`` message label.
   :param msg_txt: the ``const char*`` message text, may contain format
                   specifiers.
   :param ...: the arguments for format specifiers in ``msg_txt``.

.. c:macro:: SUNLogExtraDebugVec(logger, label, vec, msg_txt, ...)

   When extra debugging logging is enabled, this macro expands to a call to
   :c:func:`SUNLogger_QueueMsg` and :c:func:`N_VPrintFile` to log an extra
   debug message and output the vector data. Otherwise, this expands to nothing.

   :param logger: the :c:type:`SUNLogger` to handle the message.
   :param label: the ``const char*`` message label.
   :param vec: the ``N_Vector`` to print.
   :param msg_txt: the ``const char*`` message text, may contain format
                   specifiers.
   :param ...: the arguments for format specifiers in ``msg_txt``.

.. c:macro:: SUNLogExtraDebugVecIf(condition, logger, label, vec, msg_txt, ...)

   When extra debugging logging is enabled, this macro expands to a conditional
   call to :c:func:`SUNLogger_QueueMsg` and :c:func:`N_VPrintFile` to log an extra
   debug message and output the vector data. Otherwise, this expands to nothing.

   :param condition: a boolean expression that determines if the log message
                     should be queued.
   :param logger: the :c:type:`SUNLogger` to handle the message.
   :param label: the ``const char*`` message label.
   :param vec: the ``N_Vector`` to print.
   :param msg_txt: the ``const char*`` message text, may contain format
                   specifiers.
   :param ...: the arguments for format specifiers in ``msg_txt``.

.. c:macro:: SUNLogExtraDebugVecArray(logger, label, nvecs, vecs, msg_txt)

   When extra debugging logging is enabled, this macro expands to a loop calling
   :c:func:`SUNLogger_QueueMsg` and :c:func:`N_VPrintFile` for each vector in
   the vector array to log an extra debug message and output the vector data.
   Otherwise, this expands to nothing.

   :param logger: the :c:type:`SUNLogger` to handle the message.
   :param label: the ``const char*`` message label.
   :param nvecs: the ``int`` number of vectors to print.
   :param vecs: the ``N_Vector*`` (vector array) to print.
   :param msg_txt: the ``const char*`` message text, must contain a format
                   specifier for the vector array index.

   .. warning::

      The input parameter ``msg_txt`` **must** include a format specifier for
      the vector array index (of type ``int``) **only** e.g.,

      .. code-block:: C

         SUNLogExtraDebugVecArray(logger, "YS-vector-array", "YS[%d](:) =", YS, 5);


Struct Accessor Macros
----------------------

Since many SUNDIALS structs use a type-erased (i.e., `void*`) "content" pointer, 
a common idiom occurring in SUNDIALS code is extracting the content, casting it to its original
type, and then accessing the struct member of interest. To ensure readability, it is 
recommended to use locally (to the source file in question) defined macros `GET_CONTENT`
and `IMPL_MEMBER` like the following example:

.. code-block:: c

   #define GET_CONTENT(S)       ((SUNAdjointCheckpointScheme_Fixed_Content)S->content)
   #define IMPL_MEMBER(S, prop) (GET_CONTENT(S)->prop)

   SUNAdjointCheckpointScheme self;
   IMPL_MEMBER(self, current_insert_step_node)   = step_data_node;
   IMPL_MEMBER(self, step_num_of_current_insert) = step_num;
