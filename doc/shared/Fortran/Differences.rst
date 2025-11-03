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

.. _Fortran.Differences:

Notable Fortran/C usage differences
===================================

While the Fortran 2003 interface to SUNDIALS closely follows the C API, some
differences are inevitable due to the differences between Fortran and C.  In
this section, we note the most critical differences. Additionally,
:numref:`Fortran.DataTypes` discusses equivalencies of data types
in the two languages.


.. _Fortran.Differences.CreatingObjects:

Creating generic SUNDIALS objects
---------------------------------

In the C API a SUNDIALS class, such as an :c:type:`N_Vector`, is actually a pointer to
an underlying C struct. However, in the Fortran 2003 interface, the derived type
is bound to the C struct, not the pointer to the struct. For example,
``type(N_Vector)`` is bound to the C struct ``_generic_N_Vector`` not the
``N_Vector`` type. The consequence of this is that creating and declaring SUNDIALS
objects in Fortran is nuanced. This is illustrated in the code snippets below:

C code:

.. sourcecode:: c

   N_Vector x;
   x = N_VNew_Serial(N, sunctx);

Fortran code:

.. sourcecode:: Fortran

   type(N_Vector), pointer :: x
   x => FN_VNew_Serial(N, sunctx)

Note that in the Fortran declaration, the vector is a ``type(N_Vector),
pointer``, and that the pointer assignment operator is then used.


.. _Fortran.Differences.ArraysAndPointers:

Arrays and pointers
-------------------

Unlike in the C API, in the Fortran 2003 interface, arrays and pointers are
treated differently when they are return values versus arguments to a function.
Additionally, pointers which are meant to be out parameters, not arrays, in the
C API must still be declared as a rank-1 array in Fortran.  The reason for this
is partially due to the Fortran 2003 standard for C bindings, and partially due
to the tool used to generate the interfaces. Regardless, the code snippets below
illustrate the differences.

C code:

.. sourcecode:: c

   N_Vector x;
   sunrealtype* xdata;
   long int leniw, lenrw;

   /* create a new serial vector */
   x = N_VNew_Serial(N, sunctx);

   /* capturing a returned array/pointer */
   xdata = N_VGetArrayPointer(x)

   /* passing array/pointer to a function */
   N_VSetArrayPointer(xdata, x)

   /* pointers that are out-parameters */
   N_VSpace(x, &leniw, &lenrw);


Fortran code:

.. sourcecode:: Fortran

   type(N_Vector), pointer :: x
   real(c_double), pointer :: xdataptr(:)
   real(c_double)          :: xdata(N)
   integer(c_long)         :: leniw(1), lenrw(1)

   ! create a new serial vector
   x => FN_VNew_Serial(x, sunctx)

   ! capturing a returned array/pointer
   xdataptr => FN_VGetArrayPointer(x)

   ! passing array/pointer to a function
   call FN_VSetArrayPointer(xdata, x)

   ! pointers that are out-parameters
   call FN_VSpace(x, leniw, lenrw)


.. _Fortran.Differences.ProcedurePointers:

Passing procedure pointers and user data
----------------------------------------

Since functions/subroutines passed to SUNDIALS will be called from within C
code, the Fortran procedure must have the attribute ``bind(C)``. Additionally,
when providing them as arguments to a Fortran 2003 interface routine, it is
required to convert a procedure's Fortran address to C with the Fortran
intrinsic ``c_funloc``.

Typically when passing user data to a SUNDIALS function, a user may simply cast
some custom data structure as a ``void*``. When using the Fortran 2003
interfaces, the same thing can be achieved. Note, the custom data structure
*does not* have to be ``bind(C)`` since it is never accessed on the C side.

C code:

.. sourcecode:: c

   MyUserData *udata;
   void *cvode_mem;

   ierr = CVodeSetUserData(cvode_mem, udata);

Fortran code:

.. sourcecode:: Fortran

   type(MyUserData) :: udata
   type(c_ptr)      :: arkode_mem

   ierr = FARKStepSetUserData(arkode_mem, c_loc(udata))

On the other hand, Fortran users may instead choose to store problem-specific
data, e.g.  problem parameters, within modules, and thus do not need the
SUNDIALS-provided ``user_data`` pointers to pass such data back to user-supplied
functions. These users should supply the ``c_null_ptr`` input for ``user_data``
arguments to the relevant SUNDIALS functions.

.. _Fortran.Differences.OptionalParameters:

Passing ``NULL`` to optional parameters
---------------------------------------

In the SUNDIALS C API some functions have optional parameters that a caller can
pass as ``NULL``. If the optional parameter is of a type that is equivalent to a
Fortran ``type(c_ptr)`` (see :numref:`Fortran.DataTypes`),
then a Fortran user can pass the intrinsic ``c_null_ptr``. However, if the
optional parameter is of a type that is not equivalent to ``type(c_ptr)``, then
a caller must provide a Fortran pointer that is dissociated. This is
demonstrated in the code example below.

C code:

.. sourcecode:: c

   SUNLinearSolver LS;
   N_Vector x, b;

   /* SUNLinSolSolve expects a SUNMatrix or NULL as the second parameter. */
   ierr = SUNLinSolSolve(LS, NULL, x, b);

Fortran code:

.. sourcecode:: Fortran

   type(SUNLinearSolver), pointer :: LS
   type(SUNMatrix), pointer       :: A
   type(N_Vector), pointer        :: x, b

   ! Disassociate A
   A => null()

   ! SUNLinSolSolve expects a type(SUNMatrix), pointer as the second parameter.
   ! Therefore, we cannot pass a c_null_ptr, rather we pass a disassociated A.
   ierr = FSUNLinSolSolve(LS, A, x, b)

.. _Fortran.Differences.NVectorArrays:

Working with ``N_Vector`` arrays
--------------------------------

Arrays of :c:type:`N_Vector` objects are interfaced to Fortran 2003 as an opaque
``type(c_ptr)``.  As such, it is not possible to directly index an array of
:c:type:`N_Vector` objects returned by the ``N_Vector`` "VectorArray" operations, or
packages with sensitivity capabilities (CVODES and IDAS).  Instead, SUNDIALS
provides a utility function ``FN_VGetVecAtIndexVectorArray`` wrapping
:c:func:`N_VGetVecAtIndexVectorArray`. The example below demonstrates accessing
a vector in a vector array.

C code:

.. sourcecode:: c

   N_Vector x;
   N_Vector* vecs;

   /* Create an array of N_Vectors */
   vecs = N_VCloneVectorArray(count, x);

   /* Fill each array with ones */
   for (int i = 0; i < count; ++i)
     N_VConst(vecs[i], 1.0);

Fortran code:

.. sourcecode:: Fortran

   type(N_Vector), pointer :: x, xi
   type(c_ptr)             :: vecs

   ! Create an array of N_Vectors
   vecs = FN_VCloneVectorArray(count, x)

   ! Fill each array with ones
   do index = 0,count-1
     xi => FN_VGetVecAtIndexVectorArray(vecs, index)
     call FN_VConst(xi, 1.d0)
   enddo

SUNDIALS also provides the functions :c:func:`N_VSetVecAtIndexVectorArray` and
:c:func:`N_VNewVectorArray` for working with ``N_Vector`` arrays, that have
corresponding Fortran interfaces ``FN_VSetVecAtIndexVectorArray`` and
``FN_VNewVectorArray``, respectively. These functions are particularly
useful for users of the Fortran interface to the
:ref:`NVECTOR_MANYVECTOR <NVectors.ManyVector>` or
:ref:`NVECTOR_MPIMANYVECTOR <NVectors.MPIManyVector>` when creating the
subvector array. Both of these functions along with
:c:func:`N_VGetVecAtIndexVectorArray` (wrapped as
``FN_VGetVecAtIndexVectorArray``) are further described in
:numref:`NVectors.Description.utilities`.

.. _Fortran.Differences.FilePointers:

Providing file pointers
-----------------------

There are a few functions in the SUNDIALS C API which take a ``FILE*`` argument.
Since there is no portable way to convert between a Fortran file descriptor and
a C file pointer, SUNDIALS provides two utility functions for creating a
``FILE*`` and destroying it. These functions are defined in the module
``fsundials_core_mod``.

.. c:function:: SUNErrCode SUNDIALSFileOpen(const char* filename, const char* mode, FILE** fp)

   The function allocates a ``FILE*`` by calling the C function ``fopen`` with
   the provided filename and I/O mode.

   :param filename: the path to the file, that should have Fortran
      type ``character(kind=C_CHAR, len=*)``.  There are two special filenames:
      ``stdout`` and ``stderr`` -- these two filenames will result in output
      going to the standard output file and standard error file, respectively.

   :param mode: the I/O mode to use for the file.  This should have the
      Fortran type ``character(kind=C_CHAR, len=*)``.  The string begins
      with one of the following characters:

      * ``r``  to open a text file for reading
      * ``r+`` to open a text file for reading/writing
      * ``w``  to truncate a text file to zero length or create it for writing
      * ``w+`` to open a text file for reading/writing or create it if it does
        not exist
      * ``a`` to open a text file for appending, see documentation of ``fopen``
        for your system/compiler
      * ``a+`` to open a text file for reading/appending, see documentation for
        ``fopen`` for your system/compiler

   :param fp: The ``FILE*`` that will be open when the function returns.
      This should be a `type(c_ptr)` in the Fortran.

   :return: A :c:type:`SUNErrCode`

   Usage example:

   .. code-block:: Fortran

      type(c_ptr) :: fp

      ! Open up the file output.log for writing
      ierr = FSUNDIALSFileOpen("output.log", "w+", fp)

      ! The C function ARKStepPrintMem takes void* arkode_mem and FILE* fp as arguments
      call FARKStepPrintMem(arkode_mem, fp)

      ! Close the file
      ierr = FSUNDIALSFileClose(fp)

   .. versionchanged:: 7.0.0

      The function signature was updated to return a `SUNErrCode` and take a `FILE**` as the last input parameter rather then return a `FILE*`.

.. c:function:: SUNErrCode SUNDIALSFileClose(FILE** fp)

   The function deallocates a C ``FILE*`` by calling the C function ``fclose``
   with the provided pointer.

   :param fp: the C ``FILE*`` that was previously obtained from ``fopen``.
        This should have the Fortran type ``type(c_ptr)``.  Note that if either
        ``stdout`` or ``stderr`` were opened using :c:func:`SUNDIALSFileOpen()`

   :return: A :c:type:`SUNErrCode`

   .. versionchanged:: 7.0.0

      The function signature was updated to return a `SUNErrCode` and the `fp` parameter was changed from `FILE*` to `FILE**`.

