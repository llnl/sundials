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

.. _SourceCode.Rules:

Coding Conventions and Rules
============================

These rules should be followed for all new code. Unfortunately, old code might
not adhere to all of these rules.


#. Identifiers should follow our :ref:`Naming Conventions <SourceCode.Naming>`.

#. Do not use language features that are not compatible with C99, C++14,
   and MSVC v1900+ (Visual Studio 2015). Examples of such features include
   variable-length arrays. Exceptions are allowed when interfacing with a
   library which requires a newer standard.

#. All new code added to SUNDIALS should be formatted with `clang-format
   <https://clang.llvm.org/docs/ClangFormat.html>`_ for C/C++, `fprettify
   <https://github.com/fortran-lang/fprettify>`_ for Fortran, `cmake-format
   <https://cmake-format.readthedocs.io>`_ for CMake, and `black
   <https://black.readthedocs.io>`_ for Python. See :ref:`SourceCode.Style` for
   details.

#. Spaces not tabs.

#. Comments should use proper spelling and grammar.

#. Following the `Google Style Guide <https://google.github.io/styleguide/cppguide.html>`_,
   TODO comments are used to note code that is "temporary, a short-term solution,
   or good-enough but not perfect."

   A consistent TODO comment format provides an easy to search for keyword with
   details on how to get more information. TODO comments should start with ``TODO``
   followed by a unique identifier, enclosed in parentheses, for the person most
   knowledgeable about the issue and a brief description of the TODO item.
   Generally, these comments should be used sparingly and are not a substitute for
   creating an issue or bug report. When applicable, the comment should include the
   relevant issue or bug report number.

   Examples:

   .. code-block:: c

     /* TODO(DJG): Update to new API in the next major release (Issue #256) */


#. If statements and loops should always have braces even if they are one line.

#. All SUNDIALS data structures should hold onto a ``SUNContext`` object. Exceptions
   are the ``SUNLogger`` and ``SUNProfiler`` classes.

#. All SUNDIALS functions should return a ``SUNErrCode``. Many older functions
   do not do this and are exceptions to the rule for backwards compatibility.
   In addition, internal helper functions may or may-not return a ``SUNErrCode``.

#. All SUNDIALS functions, with the exception of some functions
   that do not have access to a ``SUNContext``, should begin with a call to
   ``SUNFunctionBegin()``. The argument to ``SUNFunctionBegin()`` is a ``SUNContext``
   which should come from the first object in the function parameter list that has a
   ``SUNContext``.  This macro is used for error handling and declares a
   local variable access via the macro ``SUNCTX_``.

   .. code-block:: c

      SUNErrCode N_VLinearCombination_Serial(int nvec, realtype* c, N_Vector* X, N_Vector z)
      {
         SUNFunctionBegin(X[0]->sunctx); // Correct

         int          i;
         sunindextype j, N;
         realtype*    zd=NULL;
         realtype*    xd=NULL;

         /* invalid number of vectors */
         SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

         // ...
      }

      SUNErrCode N_VLinearCombination_Serial(int nvec, realtype* c, N_Vector* X, N_Vector z)
      {
         int          i;
         sunindextype j, N;
         realtype*    zd=NULL;
         realtype*    xd=NULL;

         SUNFunctionBegin(X[0]->sunctx); // Incorrect, SUNFunctionBegin should occur as early as possible

         /* invalid number of vectors */
         SUNAssert(nvec >= 1, SUN_ERR_ARG_OUTOFRANGE);

         // ...
      }


      int CVodeGetEstLocalErrors(void *cvode_mem, N_Vector ele)
      {
         CVodeMem cv_mem;

         if (cvode_mem==NULL) {
            cvProcessError(NULL, CV_MEM_NULL, __LINE__, __func__, __FILE__, MSGCV_NO_MEM);
            return(CV_MEM_NULL);
         }

         cv_mem = (CVodeMem) cvode_mem;

         SUNFunctionBegin(cv_mem->sunctx); // Correct - this is as early as possible to call SUNFunctionBegin
         SUNFunctionBegin(ele->sunctx); // Incorrect - cvode_mem is first in the function parameter list

         // ...
      }


#. All references to ``SUNContext`` objects should be done via the ``SUNCTX_``
   macro. The only exceptions are functions in the ``SUNContext`` class.

#. All calls to SUNDIALS functions that return a ``SUNErrCode`` should have
   their return value checked with a macro from the ``SUNCheckCall`` family.
   These macros are documented in the header file ``sundials/priv/sundials_errors_impl.h``
   and ``sundials/priv/sundials_mpi_errors_impl.h``.

   .. code-block:: c

    SUNCheckCall(N_VLinearCombination(...)); // Correct

   Avoid storing the return value and then checking the stored value except when absolutely necessary.

   .. code-block:: c

    SUNErrCode err;
    err = N_VLinearCombination(...);
    SUNCheckCall(err); // Avoid except when absolutely necessary.

#. All calls to SUNDIALS functions that *do not* return a ``SUNErrCode`` should
   be followed by checking the last error stored in the ``SUNContext``.
   The exception to this rule is for internal helper functions.
   These should not be checked unless they return a ``SUNErrCode``.
   These checks are done with the ``SUNCheckLastErr`` macros.

   .. code-block:: c

    // Correct
    (void) N_VLinearSum(...); SUNCheckLastErr();

    // Incorrect - we must check for errors in N_VDotProd before calling a second function
    sunrealtype norm2 = SUNRsqrt(N_VDotProd(...)); SUNCheckLastErr();

    // Correct
    sunrealtype norm2 = N_VDotProd(...); SUNCheckLastErr();
    norm2 = SUNRsqrt(norm2);

#. Programmer errors should be checked with the ``SUNAssert`` macro, that verifies whether its
   argument evaluates to "true", and specifies an error flag otherwise. By programmer errors we
   mean, for example, illegal inputs such as mismatching dimensions or a ``NULL`` value for
   something that should not be.

   .. code-block:: c

      SUNLinearSolver SUNLinSol_Band(N_Vector y, SUNMatrix A, SUNContext sunctx)
      {
         SUNFunctionBegin(sunctx);
         SUNLinearSolver S;
         SUNLinearSolverContent_Band content;
         sunindextype MatrixRows;

         // Correct - check these with SUNAssert
         SUNAssert(SUNMatGetID(A) == SUNMATRIX_BAND, SUN_ERR_ARG_WRONGTYPE);
         SUNAssert(SUNBandMatrix_Rows(A) == SUNBandMatrix_Columns(A), SUN_ERR_ARG_DIMSMISMATCH);
         SUNAssert(y->ops->nvgetarraypointer, SUN_ERR_ARG_INCOMPATIBLE);

         // ...
      }

#. Return statements should not unnecessarily use parentheses. Prefer ``return
   x;`` to ``return(x);``. Note, however, lots of older SUNDIALS source code
   uses ``return(x);``.

#. Always use ``sunindextype`` for variables that are related to problem dimensions.
   E.g., use it for the length of a vector, or dimensions of a matrix.
   The only exception is when interfacing with a third party library requires a different
   variable type.

#. Conversely, never use ``sunindextype`` for variables that are not specifically related to
   the dimensions of a vector, matrix, etc.. E.g., if you have a variable that
   represents the number of integer "words" allocated in a workspace do not use
   ``sunindextype`` for it. Instead use the appropriate integer type (e.g., ``int64_t``) directly.
   For counters, use ``suncountertype``.

#. Do not use unsigned integer types except for ``size_t`` when the value you are storing
   is a memory size. Unsigned integer types must never be used in parts of the
   SUNDIALS API that will be interfaced to Fortran since the Fortran standard does
   not include unsigned integers.

#. Use ``sizeof(variable)`` rather than ``sizeof(type)``. E.g.,

   .. code-block:: c

      int a = 1;
      int array_length = 10;
      int* array1 = malloc(array_length * sizeof(a)); // Do this
      int* array2 = malloc(array_length * sizeof(int)); // Don't do this

#. Do not use anonymous ``enum`` s in public header files (the Python interface
   generator doesn't like it). Wrap typedef statements in SWIG guards, e.g.

   .. code-block:: c

      // Don't do this
      typedef enum {
         ARK_RELAX_BRENT,
         ARK_RELAX_NEWTON
      } ARKRelaxSolver;

      // Do this
      enum ARKRelaxSolver {
         ARK_RELAX_BRENT,
         ARK_RELAX_NEWTON
      };

      #ifndef SWIG
      typedef enum ARKRelaxSolver ARKRelaxSolver;
      #endif

#. Use the print functions, format macros, and output guidelines detailed in
   :ref:`Style.Output`.

#. Follow the logging style detailed in :ref:`Style.Logging`.
